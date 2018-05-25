#include "fully_connected_layer.h"
#include "gcop_conversions.h"
#include "load_eigen_matrix.h"
#include <Eigen/Dense>
#include <sys/stat.h>

using namespace gcop;

FullyConnectedLayer::FullyConnectedLayer(std::string variable_folder_path,
                                         std::string scope_name,
                                         std::string layer_prefix,
                                         bool use_batch_normalization,
                                         Activation activation,
                                         double batch_norm_eps)
    : use_batch_normalization_(use_batch_normalization),
      batch_norm_eps_(batch_norm_eps), activation_(activation) {
  loadParameters(variable_folder_path, layer_prefix, scope_name);
}

FullyConnectedLayer::FullyConnectedLayer() {}

cs::MX FullyConnectedLayer::transform(casadi::MX x_in) {
  cs::MX y = linearTransform(x_in, weights_, biases_);
  cs::MX yout = y;
  if (use_batch_normalization_) {
    yout = batchNorm(y, gamma_, beta_, moving_average_, moving_variance_,
                     batch_norm_eps_);
  }
  cs::MX out = activation(yout, activation_);
  return out;
}

cs::MX FullyConnectedLayer::linearTransform(cs::MX x, cs::MX weights,
                                            cs::MX biases) {
  cs::MX out = cs::MX::mtimes(weights, x);
  if (!use_batch_normalization_) {
    out = out + biases;
  }
  return out;
}

cs::MX FullyConnectedLayer::batchNorm(cs::MX x, cs::MX gamma, cs::MX beta,
                                      cs::MX moving_average,
                                      cs::MX moving_variance,
                                      double batch_norm_eps) {
  cs::MX gamma_inv_stdev =
      gamma / (cs::MX::sqrt(moving_variance + batch_norm_eps));
  cs::MX shift = (beta - gamma_inv_stdev * moving_average);
  cs::MX scaled_y = (gamma_inv_stdev)*x + shift;
  return scaled_y;
}

cs::MX FullyConnectedLayer::activation(casadi::MX x_in, Activation activation) {
  cs::MX out;
  switch (activation) {
  case Activation::none:
    out = x_in;
    break;
  case Activation::tanh:
    out = cs::MX::tanh(x_in);
    break;
  case Activation::relu:
    int rows = x_in.rows();
    int cols = x_in.columns();
    out = cs::MX::fmax(cs::MX::zeros(rows, cols), x_in);
  }
  return out;
}

void FullyConnectedLayer::loadParameters(std::string variable_folder_path,
                                         std::string layer_prefix,
                                         std::string scope_name) {
  struct stat info;
  // If folder path does not exist or is not a directory
  if (stat(variable_folder_path.c_str(), &info) != 0 ||
      !(info.st_mode & S_IFDIR)) {
    throw std::runtime_error("Cannot open folder: " + variable_folder_path);
  }
  std::string weight_file_path = addPrefixToFilePath(
      scope_name, layer_prefix, "weights", variable_folder_path);
  weights_ = loadDMFromFile(weight_file_path)
                 .T(); // Transposing weights to match shape
  if (!use_batch_normalization_) {
    std::string biases_file_path = addPrefixToFilePath(
        scope_name, layer_prefix, "biases", variable_folder_path);
    biases_ = loadDMFromFile(biases_file_path);
  } else {
    std::string beta_file_path =
        addPrefixToFilePath(scope_name, layer_prefix,
                            (scope_name + "_" + "beta"), variable_folder_path);
    std::string gamma_file_path =
        addPrefixToFilePath(scope_name, layer_prefix,
                            (scope_name + "_" + "gamma"), variable_folder_path);
    std::string moving_average_file_path = addPrefixToFilePath(
        scope_name, layer_prefix, (scope_name + "_" + "moving_mean"),
        variable_folder_path);
    std::string moving_variance_file_path = addPrefixToFilePath(
        scope_name, layer_prefix, (scope_name + "_" + "moving_variance"),
        variable_folder_path);

    beta_ = loadDMFromFile(beta_file_path);
    gamma_ = loadDMFromFile(gamma_file_path);
    moving_average_ = loadDMFromFile(moving_average_file_path);
    moving_variance_ = loadDMFromFile(moving_variance_file_path);
  }
}

std::string FullyConnectedLayer::addPrefixToFilePath(std::string scope_name,
                                                     std::string layer_prefix,
                                                     std::string variable_name,
                                                     std::string folder_path) {
  return std::string(folder_path + "/" + scope_name + "_dense_" + layer_prefix +
                     "_" + variable_name + "_0");
}

casadi::DM FullyConnectedLayer::loadDMFromFile(std::string file_path) {
  Eigen::MatrixXd mat = loadEigenMatrix(file_path);
  return conversions::convertEigenToDM(mat);
}
