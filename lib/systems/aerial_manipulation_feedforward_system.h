#ifndef AERIALMANIPULATIONFEEDFORWARDSYSTEM_H
#define AERIALMANIPULATIONFEEDFORWARDSYSTEM_H
#include "quad_casadi_system.h"

namespace gcop {

/**
 * @brief aerial manipulation feedforward model
 *
 * Combines quadrotor feedforward model with a second order model
 * for joint angles
 */
class AerialManipulationFeedforwardSystem : public CasadiSystem<> {
private:
  QuadCasadiSystem quad_system_; ///< Feedforward Quad system
  Rn<> state_manifold_;          ///< The dynamic state manifold
  Vector2d kp_ja_;               ///< Proportional Joint gains
  Vector2d kd_ja_;               ///< Derivative Joint gain
  double max_joint_velocity_;    ///< Max joint velocity
public:
  struct JointStates {
    cs::MX joint_angles;
    cs::MX joint_velocities;
    cs::MX joint_angles_desired;
  };

  /**
   * @brief Constructor
   * @param parameters The default parameter vector (thrust gain only for now)
   * @param kp_rpy The proportional gains for orientation controller
   * @param kd_rpy The derivative gains for orienation controller
   * @param kp_ja The proportional gains for joint angle controller
   * @param kd_ja The derivative gains for joint angle controller
   * @param use_code_generation Flag to specify the use of on the fly code
   * generation
   */
  AerialManipulationFeedforwardSystem(VectorXd parameters, Vector3d kp_rpy,
                                      Vector3d kd_rpy, Vector2d kp_ja,
                                      Vector2d kd_ja, double max_joint_velocity,
                                      bool use_code_generation = false);

  /**
   * @brief Aerial manipulation step function.
   *
   * Combines quad feedforward model with 2nd order joint dynamics.
   * @param t The current time
   * @param h The time step
   * @param xa The current state i.e (quad state(15) + joint state (6)).
   * @param u  The current controls i.e (quad controls (4) + joint controls(2))
   * @param p  The parameters are the quad parameters i.e thrust gain
   * @return The next state
   */
  cs::MX casadiStep(cs::MX t, cs::MX h, cs::MX xa, cs::MX u, cs::MX p);

  JointStates generateJointStates(cs::MX x);

  cs::MX jointInputs(JointStates &joint_states,
                     cs::MX joint_velocities_desired);

  cs::MX secondOrderStateUpdate(cs::MX h, JointStates &joint_states,
                                cs::MX joint_velocities_desired,
                                cs::MX joint_accelerations);

  /**
   * @brief casadiStepName
   * @return The step name
   */
  std::string casadiStepName();
};
}

#endif // AERIALMANIPULATIONFEEDFORWARDSYSTEM_H
