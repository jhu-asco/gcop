#include "aerial_manipulation_feedforward_system.h"

using namespace gcop;
using namespace casadi;

AerialManipulationFeedforwardSystem::AerialManipulationFeedforwardSystem(
    VectorXd parameters, Vector3d kp_rpy, Vector3d kd_rpy, Vector2d kp_ja,
    Vector2d kd_ja, double max_joint_velocity, bool use_code_generation)
    : CasadiSystem<>(state_manifold_, parameters, 6, 1, true, false,
                     use_code_generation),
      quad_system_(parameters, kp_rpy, kd_rpy, false), state_manifold_(21),
      kp_ja_(kp_ja), kd_ja_(kd_ja), max_joint_velocity_(max_joint_velocity) {}

AerialManipulationFeedforwardSystem::JointStates
AerialManipulationFeedforwardSystem::generateJointStates(MX x) {
  std::vector<MX> x_splits = MX::vertsplit(x, {0, 2, 4, 6});
  JointStates states;
  states.joint_angles = x_splits.at(0);
  states.joint_velocities = x_splits.at(1);
  states.joint_angles_desired = x_splits.at(2);
  return states;
}

MX AerialManipulationFeedforwardSystem::jointInputs(
    AerialManipulationFeedforwardSystem::JointStates &joint_states,
    MX joint_velocities_desired) {
  // Gains
  DM kp =
      conversions::convertEigenToDM(kp_ja_); // Proportional gains joint angles
  DM kd = conversions::convertEigenToDM(
      kd_ja_); // Derivative gains on joint velocities
  // Joint PID control
  MX e_joint_angles =
      (joint_states.joint_angles_desired - joint_states.joint_angles);
  MX e_joint_velocities =
      (joint_velocities_desired - joint_states.joint_velocities);
  MX joint_accelerations = kp * e_joint_angles + kd * e_joint_velocities;
  return joint_accelerations;
}

MX AerialManipulationFeedforwardSystem::secondOrderStateUpdate(
    MX h, AerialManipulationFeedforwardSystem::JointStates &joint_states,
    MX joint_velocities_desired, MX joint_accelerations) {
  MX joint_velocities_next =
      joint_states.joint_velocities + h * joint_accelerations;
  //  Clip joint velocities
  cs::MX constant_bnd = cs::DM::ones(2);
  joint_velocities_next =
      MX::fmax(-max_joint_velocity_ * constant_bnd, joint_velocities_next);
  joint_velocities_next =
      MX::fmin(max_joint_velocity_ * constant_bnd, joint_velocities_next);
  // clip joint velocities
  MX joint_angles_next =
      joint_states.joint_angles +
      0.5 * h * (joint_velocities_next + joint_states.joint_velocities);
  MX joint_angles_desired_next =
      joint_states.joint_angles_desired + h * joint_velocities_desired;
  return vertcat(std::vector<MX>{joint_angles_next, joint_velocities_next,
                                 joint_angles_desired_next});
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
  JointStates joint_states = generateJointStates(x_splits.at(1));
  // Quad dynamics
  MX quad_xb = quad_system_.casadiStep(t, h, quad_states, quad_controls, p);
  // Joint dynamics
  MX joint_accelerations = jointInputs(joint_states, joint_velocities_desired);
  MX joint_xb = secondOrderStateUpdate(
      h, joint_states, joint_velocities_desired, joint_accelerations);

  // Resulting state
  MX xb =
      vertcat(std::vector<MX>{quad_xb, joint_xb});
  return xb;
}

std::string AerialManipulationFeedforwardSystem::casadiStepName() {
  return "airm_ff";
}
