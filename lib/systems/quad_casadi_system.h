#ifndef QUADCASADISYSTEM_H
#define QUADCASADISYSTEM_H
#include "casadi_system.h"
#include "rn.h"

namespace gcop {
/**
 * Namespace for Eigen
 */
using namespace Eigen;
/**
 * @brief A simple feedforward model for Quad system with simple second order
 * rotational dynamics. The state of the system is given by
 * (xyz, rpy, vxyz, rpydot, rpy_desired) and control is given by (ut,
 * rpy_desired_dot)
 * The parameters for the system are given by (kt);
 */
class QuadCasadiSystem : public CasadiSystem<> {
private:
  Rn<> state_manifold_;     ///< The dynamic state manifold
  bool generate_gradients_; ///< Flag to specify whether generate gradients
  Vector3d kp_rpy_;///< Proportional gain for rpy controller
  Vector3d kd_rpy_;///< Derivative gain for rpy controller

public:
  /**
   * @brief QuadCasadiSystem Constructor
   * @param parameters The default system parameters (kt, kprpy, kdrpy)
   * @param use_code_generation Flag to decide
   * @param generate_gradients Whether to generate gradients from the step
   * function.
   */
  QuadCasadiSystem(VectorXd parameters, Vector3d kp_rpy, Vector3d kd_rpy,
                   bool use_code_generation = false,
                   bool generate_gradients = true);
  /**
   * @brief casadi_step symbolic function that integrates the system and
   * produces gradients
   * @return casadi function
   */
  cs::Function casadi_step();
  /**
   * @brief computeBodyZAxes Helper function to compute body z axis given the
   * roll pitch and yaw of the system
   * @return A function to compute body z axis given rpy of the system
   */
  cs::Function computeBodyZAxes();
};
}

#endif // QUADCASADISYSTEM_H
