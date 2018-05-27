#ifndef FULLYCONNECTEDLAYER_H
#define FULLYCONNECTEDLAYER_H
#include <casadi/casadi.hpp>
#include <gtest/gtest.h>

namespace gcop {
/**
 * @brief The namespace for casadi elements
 */
namespace cs = casadi;

/**
 * @brief The activation function used in fully connected layer
 */
enum class Activation { none, tanh, relu };

/**
 * @brief A fully connected layer similar to tensorflow model Check out
 * [link](https://www.tensorflow.org/api_docs/python/tf/contrib/layers/fully_connected)
 *
 * The fully connected layer maps an input vector into an output vector by
 * performing the following
 * functions in the order specified:
 *  - Linear transform on the input vector. Biases are not considered if using
 *    batch normalization
 *  - Batch normalization [link](https://arxiv.org/abs/1502.03167)
 *  - Activation function
 *
 */
class FullyConnectedLayer {
public:
  /**
   * @brief FullyConnectedLayer Constructor
   *
   * Loads the weights and saves parameters. The variables are saved as
   * full_name =
   * [scope_name + "_dense_" + layer_prefix + "_" + variablename +"_0"].
   *
   * For example the weights are
   * stored as [residual_dynamics_dense_1_weights_0] for scope name of
   * residual_dynamics and
   * layer prefix 1 and variable name of weights.
   *
   * @param variable_folder_path The folder to load the weights, biases and
   * batch normalization parameters
   * @param scope_name  The scope name for variables
   * @param layer_prefix  The prefix used for layers
   * @param use_batch_normalization Flag to specify whether to perform batch
   * normalization after linear transform
   * @param activation The activation to perform on the output
   * @param batch_norm_eps A small constant to add to variance to stabilize
   * batch normalization
   */
  FullyConnectedLayer(std::string variable_folder_path, std::string scope_name,
                      std::string layer_prefix,
                      bool use_batch_normalization = true,
                      Activation activation = Activation::tanh,
                      double batch_norm_eps = 0.001);
  /**
   * @brief transform input vector to output vector
   *
   * Perform linear transform, batch norm (optional) and activation
   * @param x_in Input vector can be a MX symbol or a DM
   * @return output vector
   */
  cs::MX transform(cs::MX x_in);

  /**
   * @brief linearTransform
   *
   * Multiply with weights and add biases (if batch normalization is not used)
   * @param x The input vector
   * @param weights  The weights matrix
   * @param biases The biases
   * @return The transformed matrix
   */
  cs::MX linearTransform(cs::MX x, cs::MX weights, cs::MX biases);

  /**
   * @brief batchNorm
   *
   * Perform batch normalization
   * @param x  Input vector
   * @param gamma  The scale on normalized inputs
   * @param beta  The  shift on normalized inputs
   * @param moving_average  The mean of the training data
   * @param moving_variance  The variance of the training data
   * @param batch_norm_eps A small constant to add to variance
   * @return The batch normalized input
   */
  cs::MX batchNorm(cs::MX x, cs::MX gamma, cs::MX beta, cs::MX moving_average,
                   cs::MX moving_variance, double batch_norm_eps);

  /**
   * @brief activation
   *
   * Perform activation on the input (if activation is not none)
   * @param x_in  The input vector
   * @param activation The activation function to apply as specified in
   * Activation enum
   * @return The activated inputs
   */
  cs::MX activation(cs::MX x_in, Activation activation);

  /**
   * @brief loadDMFromFile
   *
   * load casadi DataMatrix from a file. First saves it as eigen matrix
   * and then copies it to DataMatrix
   *
   * @param file_path The data file
   * @return The loaded casadi DataMatrix
   */
  cs::DM loadDMFromFile(std::string file_path);

private:
  /**
   * @brief FullyConnectedLayer
   * Private Constructor for testing private functions
   */
  FullyConnectedLayer() = default;
  /**
   * @brief loadParameters
   *
   * Load weights, biases and batch normalization parameters from a folder
   * @param variable_folder The folder where weights are present as files
   * @param layer_prefix  The prefix of the layer
   * @param scope_name The variable scope
   */
  void loadParameters(std::string variable_folder, std::string layer_prefix,
                      std::string scope_name);
  /**
   * @brief addPrefixToFilePath
   *
   * Helper function to find the path to a variable
   * @param scope_name Scope of the variable
   * @param layer_prefix  Prefix of the layer
   * @param variable_name  Name of the variable
   * @param folder_path  Path to the folder where variable is present
   * @return  full path to the variable file
   */
  std::string addPrefixToFilePath(std::string scope_name,
                                  std::string layer_prefix,
                                  std::string variable_name,
                                  std::string folder_path);
  /**
   * @brief FRIEND_TEST for testing addPrefixToFilePath
   */
  FRIEND_TEST(TestFullyConnectedLayer, addPrefixToFilePath);
  /**
   * @brief FRIEND_TEST for testing parameter loading
   */
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
