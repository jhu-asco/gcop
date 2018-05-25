#include "fully_connected_layer.h"
#include "gcop_conversions.h"
#include "load_eigen_matrix.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>

#include <chrono>

namespace gcop {

TEST(TestFullyConnectedLayer, addPrefixToFilePath) {
  FullyConnectedLayer layer;
  std::string weights_file_name = layer.addPrefixToFilePath(
      "residual_dynamics", "1", "weights", "my_folder");
  ASSERT_STREQ(weights_file_name.c_str(),
               "my_folder/residual_dynamics_dense_1_weights_0");
  std::string biases_file_name = layer.addPrefixToFilePath(
      "residual_dynamics", "final", "biases", "my_folder");
  ASSERT_STREQ(biases_file_name.c_str(),
               "my_folder/residual_dynamics_dense_final_biases_0");
}

TEST(TestFullyConnectedLayer, checkParameters) {
  std::string folder_path =
      (std::string(DATA_PATH) + "/model_vars_testing_fc_layer/");
  FullyConnectedLayer layer(folder_path, "residual_dynamics", "1", true,
                            Activation::tanh);
  ASSERT_EQ(layer.weights_.rows(), 32);
  ASSERT_EQ(layer.weights_.columns(), 36);
  ASSERT_EQ(layer.gamma_.rows(), 32);
  ASSERT_EQ(layer.beta_.rows(), 32);
  ASSERT_EQ(layer.moving_average_.rows(), 32);
  ASSERT_EQ(layer.moving_variance_.rows(), 32);
  ASSERT_FLOAT_EQ((double)layer.gamma_(0, 0), 1.0377356);
  ASSERT_FLOAT_EQ((double)layer.weights_(6, 0), 0.07734999);
  ASSERT_FLOAT_EQ((double)layer.gamma_(6, 0), 1.0190463);
}

TEST(TestFullyConnectedLayer, testTransform) {
  std::string folder_path =
      (std::string(DATA_PATH) + "/model_vars_testing_fc_layer/");
  FullyConnectedLayer layer(folder_path, "residual_dynamics", "1", true,
                            Activation::tanh);
  cs::MX x_in = cs::MX::sym("in_state", 36);
  cs::MX x_out = layer.transform(x_in);
  auto f = cs::Function("fc", {x_in}, {x_out});
  std::vector<cs::DM> args;
  args.resize(1);
  Eigen::MatrixXd inputs =
      loadEigenMatrix(folder_path + "/fc_in_0").transpose();
  Eigen::MatrixXd outputs =
      loadEigenMatrix(folder_path + "/residual_dynamics_dense_1_Tanh_0")
          .transpose();
  int N = inputs.cols();
  for (int i = 0; i < N; ++i) {
    args[0] = conversions::convertEigenToDM(inputs.col(i));
    auto out_args = f(args);
    Eigen::VectorXd cs_out = conversions::convertDMToEigen(out_args.at(0));
    Eigen::VectorXd tf_out = outputs.col(i);
    ASSERT_EQ(cs_out.rows(), tf_out.rows());
    int Nr = cs_out.rows();
    for (int j = 0; j < Nr; ++j) {
      ASSERT_NEAR(cs_out(j), tf_out(j), 1e-6);
    }
  }
}
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
