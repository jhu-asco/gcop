#include "airm_residual_network_model.h"
#include <Eigen/Dense>

using namespace gcop;
using namespace Eigen;
using namespace casadi;

AirmResidualNetworkModel::AirmResidualNetworkModel(
    VectorXd parameters, Vector3d kp_rpy, Vector3d kd_rpy, Vector2d kp_ja,
    Vector2d kd_ja, double max_joint_velocity, int n_layers,
    std::string nn_weights_folder_path, VectorXd lb, VectorXd ub,
    gcop::Activation activation, bool use_code_generation, double yaw_offset)
    : AirmResidualNetworkModel(
          parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, max_joint_velocity,
          AirmResidualNetworkModel::loadFullyConnectedLayers(
              n_layers, nn_weights_folder_path, activation),
          lb, ub, use_code_generation, yaw_offset) {}

AirmResidualNetworkModel::AirmResidualNetworkModel(
    VectorXd parameters, Vector3d kp_rpy, Vector3d kd_rpy, Vector2d kp_ja,
    Vector2d kd_ja, double max_joint_velocity,
    std::vector<FullyConnectedLayer> nn_layers, VectorXd lb, VectorXd ub,
    bool use_code_generation, double yaw_offset)
    : CasadiSystem<>(state_manifold_, parameters, 6, 1, true, false,
                     use_code_generation),
      state_manifold_(21),
      airm_system_(parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, max_joint_velocity,
                   lb, ub, false),
      quad_system_(parameters, kp_rpy, kd_rpy, lb.segment<6>(0),
                   ub.segment<6>(0), false),
      nn_layers_(nn_layers), yaw_offset_(yaw_offset) {
  std::cout << "Setting control bounds" << std::endl;
  // Set control bounds
  (this->U).lb = lb;
  (this->U).ub = ub;
  (this->U).bnd = true;
}

std::vector<FullyConnectedLayer>
AirmResidualNetworkModel::loadFullyConnectedLayers(
    int n_layers, std::string nn_weights_folder_path, Activation activation) {
  std::vector<FullyConnectedLayer> nn_layers;
  std::cout << "Loading layers..." << std::endl;
  if (n_layers < 1) {
    throw std::runtime_error(
        "The number of layers should be greater than equal to 1.\n "
        "The last layer corresponds to linear transform "
        "without any activation function/batch "
        "normalization");
  }
  for (int i = 0; i < n_layers - 1; ++i) {
    std::cout << "Loading layer: " << i << std::endl;
    nn_layers.emplace_back(nn_weights_folder_path, "residual_dynamics",
                           std::to_string(i), true, activation);
  }
  std::cout << "Loading final layer" << std::endl;
  nn_layers.emplace_back(nn_weights_folder_path, "residual_dynamics", "final",
                         false, Activation::none);
  std::cout << "Done loading layers" << std::endl;
  return nn_layers;
}

MX AirmResidualNetworkModel::propagateNetwork(MX &input) {
  MX output = input;
  for (auto &layer : nn_layers_) {
    output = layer.transform(output);
  }
  return output;
}

void AirmResidualNetworkModel::generateResidualInputs(
    QuadInputs &quad_inputs, MX &joint_accelerations, const MX &quad_states,
    const MX &joint_states, const MX &controls, const MX &p,
    const double &yaw_offset) {
  // Update Feedforward using network
  // ff inputs   + state+kt     + controls
  // 3 + 3 + 2   + 13 + 1 + 6   + 6 = 34
  MX network_vec = vertcat(std::vector<MX>{
      quad_inputs.acc, quad_inputs.rpy_ddot, joint_accelerations,
      quad_states(cs::Slice(2, 15)), p, joint_states, controls});
  MX res_vec =
      propagateNetwork(network_vec); // split into delta quad feedforward
                                     // inputs and joint accelerations
  // Update Feedforward using output from residual network
  quad_inputs.acc = quad_inputs.acc + res_vec(cs::Slice(0, 3));
  quad_inputs.rpy_ddot = quad_inputs.rpy_ddot + res_vec(cs::Slice(3, 6));
  joint_accelerations = joint_accelerations + res_vec(cs::Slice(6, 8));
}

MX AirmResidualNetworkModel::casadiStep(const MX &, const MX &h, const MX &xa,
                                        const MX &u, const MX &p) {
  // Controls [quad_controls, jad_dot]
  QuadControls quad_controls =
      quad_system_.generateControls(u(cs::Slice(0, 4)));
  const MX &joint_velocities_desired = u(cs::Slice(4, 6));
  // State
  const MX &x0 = xa(cs::Slice(0, 15));
  const MX &x1 = xa(cs::Slice(15, 21));
  QuadStates quad_states = quad_system_.generateStates(x0);
  JointStates joint_states = airm_system_.generateJointStates(x1);
  // Feedforward
  MX joint_accelerations =
      airm_system_.jointInputs(joint_states, joint_velocities_desired);
  QuadInputs quad_feedforward_inputs =
      quad_system_.computeFeedforwardInputs(quad_states, quad_controls);

  generateResidualInputs(quad_feedforward_inputs, joint_accelerations, x0, x1,
                         u, p, yaw_offset_);
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
