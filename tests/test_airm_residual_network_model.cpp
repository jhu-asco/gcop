#include "airm_residual_network_model.h"
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <memory>

using namespace gcop;

class TestAirmResidualNetworkModel : public testing::Test {
protected:
  virtual void SetUp() {
    parameters.resize(1);
    parameters << 0.18;
    kp_rpy << 15.25, 15.71, 0;
    kd_rpy << 8, 6, 3.5;
    kp_ja << 11, 11;
    kd_ja << 5, 10;
    std::string folder_path =
        (std::string(DATA_PATH) + "/tensorflow_model_vars/");
    airm_model.reset(new AirmResidualNetworkModel(
        parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, 0.7, 3, folder_path,
        Activation::tanh, true));
    airm_model->instantiateStepFunction();
  }
  template <class T> void assertVector(casadi::DM input, T expected_value) {
    std::vector<double> elems = input.get_elements();
    ASSERT_EQ(elems.size(), expected_value.size());
    ASSERT_GT(elems.size(), 0);
    for (int i = 0; i < elems.size(); ++i) {
      ASSERT_NEAR(elems[i], expected_value[i], 1e-7);
    }
  }
  std::unique_ptr<AirmResidualNetworkModel> airm_model;
  Eigen::VectorXd parameters;
  Eigen::Vector3d kp_rpy;
  Eigen::Vector3d kd_rpy;
  Eigen::Vector2d kp_ja;
  Eigen::Vector2d kd_ja;
};

TEST_F(TestAirmResidualNetworkModel, TestStep) {
  Eigen::VectorXd xa(21);
  // p,                rpy,         v,         rpydot,     rpyd,          ja,
  // jv,        jad;
  xa << 0.1, 0.1, 0.1, 0, 0, 0, 0.1, -0.2, -0.1, 0, 0, 0, 0.2, -0.2, 0.5, 0.1,
      0.2, -0.1, 0.0, 0.5, 0.5;
  Eigen::VectorXd xb(21);
  Eigen::VectorXd u(6);
  Eigen::MatrixXd A, B;
  u << 9.81 / 0.18, 0, 0, 0.1, 0, 0;
  double h = 0.02;
  int N = 3.0 / h; // tf/h
  LoopTimer overall_loop_timer;
  for (int i = 0; i < N; ++i) {
    overall_loop_timer.loop_start();
    airm_model->Step(xb, 0, xa, u, h, 0, &A, &B);
    xa = xb;
    overall_loop_timer.loop_end();
  }
  std::cout << "Copy loop average time: "
            << (airm_model->copy_loop_timer_).average_loop_period()
            << std::endl;
  std::cout << "Function loop average time: "
            << (airm_model->fun_loop_timer_).average_loop_period() << std::endl;
  std::cout << "Overall loop average time: "
            << (overall_loop_timer).average_loop_period() << std::endl;
  // rp
  ASSERT_NEAR(xb(3), 0.2, 0.1);
  ASSERT_NEAR(xb(4), -0.2, 0.1);
  // rpydot
  ASSERT_NEAR(xb(9), 0.0, 0.1);
  ASSERT_NEAR(xb(10), 0.0, 0.1);
  ASSERT_NEAR(xb(11), 0.1, 0.1);
  // rpd
  ASSERT_DOUBLE_EQ(xb(12), 0.2);
  ASSERT_DOUBLE_EQ(xb(13), -0.2);
  // ja
  ASSERT_NEAR(xb(15), 0.5, 0.1);
  ASSERT_NEAR(xb(16), 0.5, 0.1);
  // jv
  ASSERT_NEAR(xb(17), 0.0, 0.1);
  ASSERT_NEAR(xb(18), 0.0, 0.1);
  // jad
  ASSERT_NEAR(xb(19), 0.5, 0.1);
  ASSERT_NEAR(xb(20), 0.5, 0.1);
  ASSERT_EQ(A.rows(), 21);
  ASSERT_EQ(A.cols(), 21);
  ASSERT_EQ(B.rows(), 21);
  ASSERT_EQ(B.cols(), 6);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
