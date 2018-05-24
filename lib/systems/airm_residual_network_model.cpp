#include "airm_residual_network_model.h"
#include <Eigen/Dense>

using namespace gcop;
using namespace Eigen;
using namespace casadi;

AirmResidualNetworkModel::AirmResidualNetworkModel(
    VectorXd parameters, Vector3d kp_rpy, Vector3d kd_rpy, Vector2d kp_ja,
    Vector2d kd_ja, int n_layers, std::string nn_weights_folder_path,
    gcop::Activation activation, bool use_code_generation)
    : CasadiSystem<>(state_manifold_, parameters, 6, 1, true, false,
                     use_code_generation),
      airm_system_(parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, false),
      quad_system_(parameters, kp_rpy, kd_rpy, false) {
  if (n_layers <= 1) {
    throw std::runtime_error("The number of layers should be greater than 1.\n "
                             "The last layer corresponds to linear transform "
                             "without any activatin function/batch "
                             "normalization");
  }
  for (int i = 0; i < n_layers; ++i) {
    nn_layers_.emplace_back(nn_weights_folder_path, "residual_dynami",
                            std::to_string(i), true, activation);
  }
  nn_layers_.emplace_back(nn_weights_folder_path, "residual_dynami", "final",
                          false, Activation::none);
}

MX AirmResidualNetworkModel::findResidualInputs(MX input) {
  cs::MX temp = input;
  for (auto &layer : nn_layers_) {
    temp = layer.transform(temp);
  }
  return temp;
}

MX AirmResidualNetworkModel::casadiStep(MX, MX h, MX xa, MX u, MX p) {
  // Controls [quad_controls, jad_dot]
  std::vector<MX> u_splits = MX::vertsplit(u, {0, 4, 6});
  QuadControls quad_controls = quad_system_.generateControls(u_splits.at(0));
  MX joint_velocities_desired = u_splits.at(1);
  // State
  std::vector<MX> x_splits = MX::vertsplit(xa, {0, 15, 21});
  QuadStates quad_states = quad_system_.generateStates(x_splits.at(0));
  JointStates joint_states = airm_system_.generateJointStates(x_splits.at(1));
  // Feedforward
  MX joint_accelerations =
      airm_system_.jointInputs(joint_states, joint_velocities_desired);
  QuadInputs quad_feedforward_inputs =
      quad_system_.computeFeedforwardInputs(quad_states, quad_controls, p);
  // Update Feedforward using network
  // ff inputs   + state+kt     + controls
  // 3 + 3 + 2   + 15 + 1 + 6   + 6 = 36
  MX network_input = vertcat(std::vector<MX>{
      quad_feedforward_inputs.acc, quad_feedforward_inputs.rpy_ddot,
      joint_accelerations, x_splits.at(0), p, x_splits.at(1), u});
  MX res_feedforward_input =
      findResidualInputs(network_input); // split into delta quad feedforward
                                         // inputs and joint accelerations
  // Update Feedforward using output from residual network
  std::vector<MX> res_input_split =
      MX::vertsplit(res_feedforward_input, {0, 3, 6, 8});
  quad_feedforward_inputs.acc =
      quad_feedforward_inputs.acc + res_input_split.at(0);
  quad_feedforward_inputs.rpy_ddot =
      quad_feedforward_inputs.rpy_ddot + res_input_split.at(1);
  joint_accelerations = joint_accelerations + res_input_split.at(2);
  // Integrate
  MX quad_xb = quad_system_.secondOrderStateUpdate(
      h, quad_states, quad_controls, quad_feedforward_inputs);
  MX joint_xb = airm_system_.secondOrderStateUpdate(
      h, joint_states, joint_velocities_desired, joint_accelerations);
  return vertcat(std::vector<MX>{quad_xb, joint_xb});
}

std::string AirmResidualNetworkModel::casadiStepName() {
  return "airm_res_step";
}
