#ifndef FULLYCONNECTEDLAYER_H
#define FULLYCONNECTEDLAYER_H
#include <casadi/casadi.hpp>
#include <gtest/gtest.h>

namespace gcop {
namespace cs = casadi;

enum class Activation { none, tanh, relu };

class FullyConnectedLayer {
public:
  FullyConnectedLayer(std::string variable_folder_path, std::string scope_name,
                      std::string layer_prefix, bool use_batch_normalization,
                      Activation activation, double batch_norm_eps = 0.001);

  FullyConnectedLayer();

  cs::MX transform(cs::MX x_in);

  cs::MX linearTransform(cs::MX x, cs::DM weights, cs::DM biases);

  cs::MX batchNorm(cs::MX x, cs::DM gamma, cs::DM beta, cs::DM moving_average,
                   cs::DM moving_variance, double batch_norm_eps);

  cs::MX activation(cs::MX x_in, Activation activation);

  cs::DM loadDMFromFile(std::string file_path);

private:
  void loadParameters(std::string variable_folder, std::string layer_prefix,
                      std::string scope_name);
  std::string addPrefixToFilePath(std::string scope_name,
                                  std::string layer_prefix,
                                  std::string variable_name,
                                  std::string folder_path);
  FRIEND_TEST(TestFullyConnectedLayer, addPrefixToFilePath);
  FRIEND_TEST(TestFullyConnectedLayer, checkParameters);

private:
  bool use_batch_normalization_; ///< Whether to use batch normalization or not
  double batch_norm_eps_;        ///< Small constant for inverting variance
  Activation activation_;        ///< Type of activation layer
  cs::DM weights_; ///< Weight matrix for transforming inputs (Transpose of the
                   /// Tensorflow matrix i.e W * x is used instead of x^T *W^T)
  cs::DM
      biases_; ///< Bias matrix used only if batch normalization is not present
  cs::DM gamma_;           ///< Scale for batch normalization
  cs::DM beta_;            ///< Bias for batch normalization
  cs::DM moving_average_;  ///< Moving average from batch normalization
  cs::DM moving_variance_; ///< Moving variance from batch normalization
};
}

#endif // FULLYCONNECTEDLAYER_H
