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
  Rn<> state_manifold_;          ///< The dynamic state manifold
  Vector2d kp_ja_;               ///< Proportional Joint gains
  Vector2d kd_ja_;               ///< Derivative Joint gain
  QuadCasadiSystem quad_system_; ///< Feedforward Quad system
public:
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
                                      Vector2d kd_ja, bool use_code_generation);

  /**
   * @brief Aerial manipulation step function.
   *
   * Combines quad feedforward model with 2nd order joint dynamics.
   * @param t The current time
   * @param h The time step
   * @param xa The current state i.e (quad state(15) + joint state (6)).
   *           The joint state consists of joint angles, velocities and desired
   *           joint angles
   * @param u  The current controls i.e (quad controls (4) + joint controls(2))
   *           The joint controls are the desired joint velocities
   * @param p  The parameters are the quad parameters i.e thrust gain
   * @return The next state
   */
  cs::MX casadi_step(cs::MX t, cs::MX h, cs::MX xa, cs::MX u, cs::MX p);
  /**
   * @brief casadi_step_name
   * @return The step name
   */
  std::string casadi_step_name();
};
}

#endif // AERIALMANIPULATIONFEEDFORWARDSYSTEM_H
