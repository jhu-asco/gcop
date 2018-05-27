#include "load_eigen_matrix.h"
#include <Eigen/Dense>
#include <gtest/gtest.h>

TEST(LoadEigen, UnknownFile) {
  ASSERT_THROW(gcop::loadEigenMatrix("unknown_file"), std::runtime_error);
}

TEST(LoadEigen, NegativeCols) {
  // Data path defined in CMakelists
  ASSERT_THROW(gcop::loadEigenMatrix(std::string(DATA_PATH) +
                                     "/test_model_vars/kd_rpy_negative_col"),
               std::runtime_error);
}

TEST(LoadEigen, WrongData) {
  ASSERT_THROW(gcop::loadEigenMatrix(std::string(DATA_PATH) +
                                     "/test_model_vars/kd_rpy_wrong_data"),
               std::runtime_error);
}

TEST(LoadEigen, VectorXd) {
  Eigen::VectorXd rpy_gains =
      gcop::loadEigenMatrix(std::string(DATA_PATH) + "/test_model_vars/kd_rpy");
  ASSERT_EQ(rpy_gains.rows(), 3);
  ASSERT_EQ(rpy_gains.cols(), 1);
  for (int i = 0; i < 3; ++i) {
    ASSERT_DOUBLE_EQ(rpy_gains[i], i + 1);
  }
}

TEST(LoadEigen, MatrixXd) {
  Eigen::MatrixXd weights = gcop::loadEigenMatrix(
      std::string(DATA_PATH) + "/test_model_vars/dense_0_weights_0");
  ASSERT_EQ(weights.rows(), 36);
  ASSERT_EQ(weights.cols(), 32);
  ASSERT_FLOAT_EQ(weights(0, 0), -0.1492475);
  ASSERT_FLOAT_EQ(weights(16, 0), 0.2718344);
  ASSERT_FLOAT_EQ(weights(0, 16), 0.23243913);
  ASSERT_FLOAT_EQ(weights(8, 12), 0.018403042);
  ASSERT_FLOAT_EQ(weights(12, 8), 0.46349743); // Verified from python
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
