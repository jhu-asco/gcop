#include "quad_casadi_system.h"

using namespace gcop;
using namespace Eigen;
using namespace casadi;

QuadCasadiSystem::QuadCasadiSystem(VectorXd parameters,
                                   bool use_code_generation,
                                   bool generate_gradients)
    : CasadiSystem<>(state_manifold_, parameters, 4, 7, use_code_generation),
      state_manifold_(15), generate_gradients_(generate_gradients) {
  // States: p, rpy, v, rpydot, rpyd
  // Controls: ut, rpyd (4)
  // Parameters kt, kp_rpy, kd_rpy
  casadi_assert(parameters.size() == 7, "Size of parameters should be 7");
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

Function QuadCasadiSystem::casadi_step() {
  MX t = MX::sym("t", 1);   // Current time
  MX dt = MX::sym("dt", 1); // time step

  // Controls
  MX ut = MX::sym("ut", 1);                      // Thrust command
  MX rpy_dot_desired = MX::sym("rpydot_des", 3); // Desired rpy dot
  MX u = vertcat(ut, rpy_dot_desired);

  // State
  MX p_i = MX::sym("p_i", 3);                     // initial position
  MX rpy_i = MX::sym("rpy_i", 3);                 // initial rpy
  MX v_i = MX::sym("v_i", 3);                     // initial velocity
  MX rpy_dot_i = MX::sym("rpy_dot_i", 3);         // initial rpydot
  MX rpy_desired_i = MX::sym("rpy_desired_i", 3); // initial rpy_desired

  MX X =
      vertcat(std::vector<MX>{p_i, rpy_i, v_i, rpy_dot_i, rpy_desired_i}); // 15

  // Internal params
  MX kt = MX::sym("kt", 1);         // thrust gain
  MX kp_rpy = MX::sym("kp_rpy", 3); // Proportional gains rpy
  MX kd_rpy = MX::sym("kd_rpy", 3); // Derivative gains on rpydot
  MX p = vertcat(kt, kp_rpy, kd_rpy);
  MX g = MX({0, 0, -9.81});

  // Constants
  Function body_z_axis_comp = computeBodyZAxes();
  MX z_axis = body_z_axis_comp(std::vector<MX>{rpy_i}).at(0);
  MX acc = kt * u(0) * z_axis + g;
  // Translational dynamics
  MX v_next = v_i + dt * acc;
  MX p_next = p_i + 0.5 * dt * (v_next + v_i);
  // Rotational dynamics
  MX e_rpy = (rpy_desired_i - rpy_i);
  MX e_rpy_dot = (rpy_dot_desired - rpy_dot_i);
  MX rpyddot = kp_rpy * e_rpy + kd_rpy * e_rpy_dot;
  MX rpy_dot_next = rpy_dot_i + dt * rpyddot;
  MX rpy_next = rpy_i + 0.5 * dt * (rpy_dot_next + rpy_dot_i);
  MX rpy_desired_next = rpy_desired_i + dt * rpy_dot_desired;
  // Resulting state
  MX Xn = vertcat(std::vector<MX>{p_next, rpy_next, v_next, rpy_dot_next,
                                  rpy_desired_next});

  if (generate_gradients_) {
    MX A = MX::jacobian(Xn, X);
    MX B = MX::jacobian(Xn, u);
    return Function("quad_step", {t, dt, X, u, p}, {Xn, A, B});
  }

  return Function("quad_step", {t, dt, X, u, p}, {Xn});
}
