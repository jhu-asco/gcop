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
public:
  struct States {
    cs::MX p;
    cs::MX rpy;
    cs::MX v;
    cs::MX rpy_dot;
    cs::MX rpy_desired;
  };

  struct Controls {
    cs::MX ut;
    cs::MX rpy_dot_desired;
  };
  struct FeedforwardInputs {
    cs::MX acc;
    cs::MX rpy_ddot;
  };

  /**
   * @brief QuadCasadiSystem Constructor
   * @param parameters The default system parameters (kt, kprpy, kdrpy)
   * @param use_code_generation Flag to decide
   * function.
   */
  QuadCasadiSystem(VectorXd parameters, Vector3d kp_rpy, Vector3d kd_rpy,
                   bool use_code_generation = false);

  FeedforwardInputs computeFeedforwardInputs(States &x_splits,
                                             Controls &u_splits);

  cs::MX secondOrderStateUpdate(const cs::MX &h, const States &x_splits,
                                const Controls &u_splits,
                                const FeedforwardInputs &feedforward_inputs);

  States generateStates(const cs::MX &x);

  Controls generateControls(const cs::MX &u);

  /**
   * @brief casadi_step
   * @param t Current time
   * @param h Time step
   * @param xa Current state
   * @param u Current control
   * @param p Parameter
   * @return Next state
   */
  cs::MX casadiStep(const cs::MX &, const cs::MX &h, const cs::MX &xa,
                    const cs::MX &u, const cs::MX &p);

  /**
  * @brief The name the step function
  *
  * @return name of step function
  */
  std::string casadiStepName();
  /**
   * @brief computeBodyZAxes Helper function to compute body z axis given the
   * roll pitch and yaw of the system
   * @return A function to compute body z axis given rpy of the system
   */
  cs::MX computeBodyZAxes(const cs::MX &rpy);

private:
  Rn<> state_manifold_; ///< The dynamic state manifold
  Vector3d kp_rpy_;     ///< Proportional gain for rpy controller
  Vector3d kd_rpy_;     ///< Derivative gain for rpy controller
};
}

#endif // QUADCASADISYSTEM_H
