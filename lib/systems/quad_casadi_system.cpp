#include "quad_casadi_system.h"

using namespace gcop;
using namespace Eigen;
using namespace casadi;

QuadCasadiSystem::QuadCasadiSystem(VectorXd parameters, Vector3d kp_rpy,
                                   Vector3d kd_rpy, bool use_code_generation)
    : CasadiSystem<>(state_manifold_, parameters, 4, 1, true, false,
                     use_code_generation),
      state_manifold_(15), kp_rpy_(kp_rpy), kd_rpy_(kd_rpy) {
  // States: p, rpy, v, rpydot, rpyd
  // Controls: ut, rpyd (4)
  casadi_assert(parameters.size() == 1, "Size of parameters should be 1");
}

Function QuadCasadiSystem::computeBodyZAxes() {
  MX rpy = MX::sym("rpy", 3);
  MX c_rpy = MX::cos(rpy);
  MX s_rpy = MX::sin(rpy);
  MX temp = c_rpy(0) * s_rpy(1);
  MX z_axis =
      vertcat(temp * c_rpy(2) + s_rpy(0) * s_rpy(2),
              temp * s_rpy(2) - s_rpy(0) * c_rpy(2), c_rpy(0) * c_rpy(1));
  return Function("compute_body_z_axes", {rpy}, {z_axis});
}

cs::MX QuadCasadiSystem::casadiStep(MX, MX h, MX xa, MX u, MX p) {

  // Controls
  // Assuming offset is [o1, o2, o3,..., on] then n-1 vectors
  // are generated of sizes [oi-o(i-1)] starting from oi.
  std::vector<MX> u_splits = MX::vertsplit(u, {0, 1, 4});
  MX ut = u_splits.at(0);              // Thrust command
  MX rpy_dot_desired = u_splits.at(1); // Desired rpy dot

  // State
  std::vector<MX> x_splits = MX::vertsplit(xa, {0, 3, 6, 9, 12, 15});
  MX p_i = x_splits.at(0);           // initial position
  MX rpy_i = x_splits.at(1);         // initial rpy
  MX v_i = x_splits.at(2);           // initial velocity
  MX rpy_dot_i = x_splits.at(3);     // initial rpydot
  MX rpy_desired_i = x_splits.at(4); // initial rpy_desired

  // Internal params
  MX kt = p; // thrust gain
  DM kp_rpy = eigen_casadi_conversions::convertEigenToDM(
      kp_rpy_); // Proportional gains rpy
  DM kd_rpy = eigen_casadi_conversions::convertEigenToDM(
      kd_rpy_); // Derivative gains on rpydot
  MX g = MX({0, 0, -9.81});

  // Constants
  Function body_z_axis_comp = computeBodyZAxes();
  MX z_axis = body_z_axis_comp(std::vector<MX>{rpy_i}).at(0);
  MX acc = kt * ut * z_axis + g;
  // Translational dynamics
  MX v_next = v_i + h * acc;
  MX p_next = p_i + 0.5 * h * (v_next + v_i);
  // Rotational dynamics
  MX e_rpy = (rpy_desired_i - rpy_i);
  MX e_rpy_dot = (rpy_dot_desired - rpy_dot_i);
  MX rpyddot = kp_rpy * e_rpy + kd_rpy * e_rpy_dot;
  MX rpy_dot_next = rpy_dot_i + h * rpyddot;
  MX rpy_next = rpy_i + 0.5 * h * (rpy_dot_next + rpy_dot_i);
  MX rpy_desired_next = rpy_desired_i + h * rpy_dot_desired;
  // Resulting state
  MX xb = vertcat(std::vector<MX>{p_next, rpy_next, v_next, rpy_dot_next,
                                  rpy_desired_next});
  return xb;
}

std::string QuadCasadiSystem::casadiStepName() { return "quad_step"; }
