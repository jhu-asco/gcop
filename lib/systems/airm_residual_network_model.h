#ifndef AIRMRESIDUALNETWORKMODEL_H
#define AIRMRESIDUALNETWORKMODEL_H
#include "aerial_manipulation_feedforward_system.h"
#include "fully_connected_layer.h"

namespace gcop {
class AirmResidualNetworkModel : public CasadiSystem<> {
private:
  /**
   * @brief The dynamic state manifold
   */
  Rn<> state_manifold_;
  /**
   * @brief Feedforward Aerial manipulation system
   */
  AerialManipulationFeedforwardSystem airm_system_;
  /**
   * @brief Internal quad system to find quad feedforward inputs
   */
  QuadCasadiSystem quad_system_;
  /**
   * @brief nn_layers_ Fully connected layers
   */
  std::vector<FullyConnectedLayer> nn_layers_;

  double yaw_offset_;
  /**
   * Struct that holds joint angles, velocities, desired angles.
   */
  using JointStates = AerialManipulationFeedforwardSystem::JointStates;
  /**
   * Struct holding quadrotor state
   */
  using QuadStates = QuadCasadiSystem::States;
  /**
   * Struct holding thrust, rpydot as commands
   */
  using QuadControls = QuadCasadiSystem::Controls;
  /**
   * Struct holding feedforward inputs i.e acc, rpyddot
   */
  using QuadInputs = QuadCasadiSystem::FeedforwardInputs;

public:
  /**
   * @brief AirmResidualNetworkModel Constructor
   * @param parameters The default thrust gain
   * @param kp_rpy Proportional gains for rpy controller
   * @param kd_rpy  Derivative gains for rpy controller
   * @param kp_ja Proportional gains for joint controller
   * @param kd_ja Derivative gains for joint controller
   * @param max_joint_velocity Max joint velocity
   * @param n_layers Number of neural network layers
   * @param nn_weights_folder_path The folder with contains nn weights, biases
   * @param activation The activation layer to use
   * @param use_code_generation To use codegeneration or not when instantiating
   * step function
   */
  AirmResidualNetworkModel(VectorXd parameters, Vector3d kp_rpy,
                           Vector3d kd_rpy, Vector2d kp_ja, Vector2d kd_ja,
                           double max_joint_velocity, int n_layers,
                           std::string nn_weights_folder_path, VectorXd lb,
                           VectorXd ub,
                           Activation activation = Activation::tanh,
                           bool use_code_generation = false,
                           double yaw_offset = 0);

  /**
   * @brief propagate the input through a series of fully connected layers
   * @param input Input vector to the network
   * @return output of the several layers
   */
  casadi::MX propagateNetwork(casadi::MX &input);

  /**
   * @brief casadiStep Step function
   * @param h The timestep
   * @param xa Starting state
   * @param u Control
   * @param p Parameter i.e thrust gain kt
   * @return  The next state
   */
  casadi::MX casadiStep(const casadi::MX &, const casadi::MX &h,
                        const casadi::MX &xa, const casadi::MX &u,
                        const casadi::MX &p);
  /**
   * @brief casadiStepName
   * @return The name of the step function
   */
  std::string casadiStepName();
  /**
   * @brief generateResidualInputs
   *
   * Updates quad inputs, joint accelerations by using a fully
   * connected network.
   * @param quad_inputs The feedforward quadrotor inputs
   * @param joint_accelerations The joint accelerations
   * @param quad_states The state of the quadrotor system
   * @param joint_states The state of the joint
   * @param controls The controls to aerial manipulation system
   * @param p
   */
  void generateResidualInputs(QuadInputs &quad_inputs,
                              casadi::MX &joint_accelerations,
                              const casadi::MX &quad_states,
                              const casadi::MX &joint_states,
                              const QuadControls &quad_controls,
                              const casadi::MX &joint_desired_velocities,
                              const casadi::MX &p,
                              const double &yaw_offset = 0);
};
}

#endif // AIRMRESIDUALNETWORKMODEL_H
