#include "quad_casadi_system.h"

using namespace gcop;
using namespace Eigen;
using namespace casadi;

QuadCasadiSystem::QuadCasadiSystem(VectorXd parameters, Vector3d kp_rpy,
                                   Vector3d kd_rpy, VectorXd lb, VectorXd ub,
                                   bool use_code_generation)
    : CasadiSystem<>(state_manifold_, parameters, 4, 1, true, false,
                     use_code_generation),
      state_manifold_(15), kp_rpy_(kp_rpy), kd_rpy_(kd_rpy) {
  // States: p, rpy, v, rpydot, rpyd
  // Controls: ut, rpyd (4)
  casadi_assert(parameters.size() == 1, "Size of parameters should be 1");
  // Set control bounds
  (this->U).lb = lb;
  (this->U).ub = ub;
  (this->U).bnd = true;
}

QuadCasadiSystem::FeedforwardInputs
QuadCasadiSystem::computeFeedforwardInputs(States &x_splits,
                                           Controls &u_splits) {
  // Constants
  DM kp_rpy = conversions::convertEigenToDM(kp_rpy_); // Proportional gains rpy
  DM kd_rpy =
      conversions::convertEigenToDM(kd_rpy_); // Derivative gains on rpydot
  MX g = MX({0, 0, -9.81});
  // Output
  FeedforwardInputs inputs;
  // Acc input
  MX z_axis = computeBodyZAxes(x_splits.rpy);
  // inputs.acc = kt * u_splits.ut * z_axis + g;
  inputs.acc = 9.81 * u_splits.ut * z_axis + g;
  // rpyddot input
  MX e_rpy = (x_splits.rpy_desired - x_splits.rpy);
  MX e_rpy_dot = (u_splits.rpy_dot_desired - x_splits.rpy_dot);
  inputs.rpy_ddot = kp_rpy * e_rpy + kd_rpy * e_rpy_dot;
  return inputs;
}

MX QuadCasadiSystem::secondOrderStateUpdate(
    const MX &h, const States &x_splits, const Controls &u_splits,
    const FeedforwardInputs &feedforward_inputs) {
  // Translational dynamics
  MX v_next = x_splits.v + h * feedforward_inputs.acc;
  MX p_next = x_splits.p + 0.5 * h * (v_next + x_splits.v);
  // Rotational dynamics
  MX rpy_dot_next = x_splits.rpy_dot + h * feedforward_inputs.rpy_ddot;
  MX rpy_next = x_splits.rpy + 0.5 * h * (rpy_dot_next + x_splits.rpy_dot);
  MX rpy_desired_next = x_splits.rpy_desired + h * u_splits.rpy_dot_desired;

  return vertcat(std::vector<MX>{p_next, rpy_next, v_next, rpy_dot_next,
                                 rpy_desired_next});
}

QuadCasadiSystem::States QuadCasadiSystem::generateStates(const MX &x) {
  // State
  std::vector<MX> x_splits = MX::vertsplit(x, {0, 3, 6, 9, 12, 15});
  States state;
  state.p = x_splits.at(0);           // initial position
  state.rpy = x_splits.at(1);         // initial rpy
  state.v = x_splits.at(2);           // initial velocity
  state.rpy_dot = x_splits.at(3);     // initial rpydot
  state.rpy_desired = x_splits.at(4); // initial rpy_desired
  return state;
}

QuadCasadiSystem::Controls QuadCasadiSystem::generateControls(const MX &u) {
  // Controls
  // Assuming offset is [o1, o2, o3,..., on] then n-1 vectors
  // are generated of sizes [oi-o(i-1)] starting from oi.
  std::vector<MX> u_splits = MX::vertsplit(u, {0, 1, 4});
  Controls controls;
  controls.ut = u_splits.at(0);              // Thrust command
  controls.rpy_dot_desired = u_splits.at(1); // Desired rpy dot
  return controls;
}

MX QuadCasadiSystem::computeBodyZAxes(const MX &rpy) {
  MX c_rpy = MX::cos(rpy);
  MX s_rpy = MX::sin(rpy);
  MX temp = c_rpy(0) * s_rpy(1);
  MX z_axis =
      vertcat(temp * c_rpy(2) + s_rpy(0) * s_rpy(2),
              temp * s_rpy(2) - s_rpy(0) * c_rpy(2), c_rpy(0) * c_rpy(1));
  return z_axis;
}

cs::MX QuadCasadiSystem::casadiStep(const MX &, const MX &h, const MX &xa,
                                    const MX &u, const MX &p) {
  States x_splits = generateStates(xa);
  Controls u_splits = generateControls(u);
  FeedforwardInputs inputs = computeFeedforwardInputs(x_splits, u_splits);
  cs::MX xb = secondOrderStateUpdate(h, x_splits, u_splits, inputs);
  return xb;
}

std::string QuadCasadiSystem::casadiStepName() { return "quad_step"; }
