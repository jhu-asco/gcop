#include "aerial_manipulation_feedforward_system.h"

using namespace gcop;
using namespace casadi;

AerialManipulationFeedforwardSystem::AerialManipulationFeedforwardSystem(
    VectorXd parameters, Vector3d kp_rpy, Vector3d kd_rpy, Vector2d kp_ja,
    Vector2d kd_ja, bool use_code_generation)
    : CasadiSystem<>(state_manifold_, parameters, 6, 1, true, false,
                     use_code_generation),
      quad_system_(parameters, kp_rpy, kd_rpy, false), state_manifold_(21),
      kp_ja_(kp_ja), kd_ja_(kd_ja) {}

MX AerialManipulationFeedforwardSystem::joint_step(MX h, MX xa, MX u) {
  std::vector<MX> x_splits = MX::vertsplit(xa, {0, 2, 4, 6});
  MX joint_angles = x_splits.at(0);
  MX joint_velocities = x_splits.at(1);
  MX joint_angles_desired = x_splits.at(2);
  // Controls
  MX joint_velocities_desired = u;
  // Gains
  DM kp = eigen_casadi_conversions::convertEigenToDM(
      kp_ja_); // Proportional gains joint angles
  DM kd = eigen_casadi_conversions::convertEigenToDM(
      kd_ja_); // Derivative gains on joint velocities
  // Joint dynamics
  MX e_joint_angles = (joint_angles_desired - joint_angles);
  MX e_joint_velocities = (joint_velocities_desired - joint_velocities);
  MX joint_accelerations = kp * e_joint_angles + kd * e_joint_velocities;
  MX joint_velocities_next = joint_velocities + h * joint_accelerations;
  MX joint_angles_next =
      joint_angles + 0.5 * h * (joint_velocities_next + joint_velocities);
  MX joint_angles_desired_next =
      joint_angles_desired + h * joint_velocities_desired;
  MX xb = vertcat(std::vector<MX>{joint_angles_next, joint_velocities_next,
                              joint_angles_desired_next});
  return xb;
}

MX AerialManipulationFeedforwardSystem::casadiStep(MX t, MX h, MX xa, MX u,
                                                   MX p) {
  // Controls [quad_controls, jad_dot]
  std::vector<MX> u_splits = MX::vertsplit(u, {0, 4, 6});
  MX quad_controls = u_splits.at(0);
  MX joint_velocities_desired = u_splits.at(1);
  // State
  std::vector<MX> x_splits = MX::vertsplit(xa, {0, 15, 21});
  MX quad_states = x_splits.at(0);
  MX joint_states = x_splits.at(1);
  // Quad dynamics
  MX quad_xb = quad_system_.casadiStep(t, h, quad_states, quad_controls, p);
  // Joint dynamics
  MX joint_xb = joint_step(h, joint_states, joint_velocities_desired);

  // Resulting state
  MX xb =
      vertcat(std::vector<MX>{quad_xb, joint_xb});
  return xb;
}

std::string AerialManipulationFeedforwardSystem::casadiStepName() {
  return "airm_ff";
}
