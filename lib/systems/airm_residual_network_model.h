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
  QuadCasadiSystem quad_system_;
  std::vector<FullyConnectedLayer> nn_layers_;
  using JointStates = AerialManipulationFeedforwardSystem::JointStates;
  using QuadStates = QuadCasadiSystem::States;
  using QuadControls = QuadCasadiSystem::Controls;
  using QuadInputs = QuadCasadiSystem::FeedforwardInputs;

public:
  AirmResidualNetworkModel(VectorXd parameters, Vector3d kp_rpy,
                           Vector3d kd_rpy, Vector2d kp_ja, Vector2d kd_ja,
                           int n_layers, std::string nn_weights_folder_path,
                           Activation activation,
                           bool use_code_generation = false);

  cs::MX findResidualInputs(cs::MX input);

  cs::MX casadiStep(casadi::MX, casadi::MX h, casadi::MX xa, casadi::MX u,
                    casadi::MX p);
  std::string casadiStepName();
};
}

#endif // AIRMRESIDUALNETWORKMODEL_H
