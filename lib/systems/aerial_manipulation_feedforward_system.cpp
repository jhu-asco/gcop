#include "aerial_manipulation_feedforward_system.h"

using namespace gcop;
using namespace casadi;

AerialManipulationFeedforwardSystem::AerialManipulationFeedforwardSystem(
    VectorXd parameters, Vector3d kp_rpy, Vector3d kd_rpy, Vector2d kp_ja,
    Vector2d kd_ja, double max_joint_velocity, VectorXd lb, VectorXd ub,
    bool use_code_generation)
    : CasadiSystem<>(state_manifold_, parameters, 6, 1, true, false,
                     use_code_generation),
      quad_system_(parameters, kp_rpy, kd_rpy, lb.segment<4>(0),
                   ub.segment<4>(0), false),
      state_manifold_(21), kp_ja_(kp_ja), kd_ja_(kd_ja),
      max_joint_velocity_(max_joint_velocity) {
  // Set control bounds
  (this->U).lb = lb;
  (this->U).ub = ub;
  (this->U).bnd = true;
}

AerialManipulationFeedforwardSystem::JointStates
AerialManipulationFeedforwardSystem::generateJointStates(const MX &x) {
  JointStates states;
  states.joint_angles = x(cs::Slice(0, 2));
  states.joint_velocities = x(cs::Slice(2, 4));
  states.joint_angles_desired = x(cs::Slice(4, 6));
  return states;
}

MX AerialManipulationFeedforwardSystem::jointInputs(
    const JointStates &joint_states, const MX &joint_velocities_desired) {
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
    const MX &h, const JointStates &joint_states,
    const MX &joint_velocities_desired, const MX &joint_accelerations) {
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

MX AerialManipulationFeedforwardSystem::casadiStep(const MX &t, const MX &h,
                                                   const MX &xa, const MX &u,
                                                   const MX &p) {
  // Controls [quad_controls, jad_dot]
  MX quad_controls = u(cs::Slice(0, 4));
  MX joint_velocities_desired = u(cs::Slice(4, 6));
  // State
  MX quad_states = xa(cs::Slice(0, 15));
  JointStates joint_states = generateJointStates(xa(cs::Slice(15, 21)));
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
