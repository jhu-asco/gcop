#include "aerial_manipulation_feedforward_system.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <memory>

#include <chrono>

using namespace gcop;

class TestAerialManipulationFeedforwardSystem : public testing::Test {
protected:
  virtual void SetUp() {
    parameters.resize(1);
    parameters << 0.16;
    kp_rpy<< 10, 10, 0;
    kd_rpy<< 5, 5, 2;
    kp_ja<<10, 15;
    kd_ja<<5, 10;
    airm_system.reset(new AerialManipulationFeedforwardSystem(
        parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, 0.7));
    airm_system->instantiateStepFunction();
  }
  template <class T> void assertVector(casadi::DM input, T expected_value) {
    std::vector<double> elems = input.get_elements();
    ASSERT_EQ(elems.size(), expected_value.size());
    ASSERT_GT(elems.size(), 0);
    for (int i = 0; i < elems.size(); ++i) {
      ASSERT_NEAR(elems[i], expected_value[i], 1e-7);
    }
  }
  std::unique_ptr<AerialManipulationFeedforwardSystem> airm_system;
  Eigen::VectorXd parameters;
  Eigen::Vector3d kp_rpy;
  Eigen::Vector3d kd_rpy;
  Eigen::Vector2d kp_ja;
  Eigen::Vector2d kd_ja;
};

TEST_F(TestAerialManipulationFeedforwardSystem, Initialize) {}

TEST_F(TestAerialManipulationFeedforwardSystem, InitializeCodeGen) {
  airm_system.reset(new AerialManipulationFeedforwardSystem(
      parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, true));
  airm_system->instantiateStepFunction();
}

TEST_F(TestAerialManipulationFeedforwardSystem, TestStep) {
  // Use code generation
  airm_system.reset(new AerialManipulationFeedforwardSystem(
      parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, 0.7, true));
  airm_system->instantiateStepFunction();
  Eigen::VectorXd xa(21);
  // p,                rpy,         v,         rpydot,     rpyd,          ja,       jv,        jad;
  xa << 0.1, 0.1, 0.1, 0, 0, 0, 0.1, -0.2, -0.1, 0, 0, 0, 0.2, -0.2, 0.5, 0.1, 0.2, -0.1, 0.0, 0.5, 0.5;
  Eigen::VectorXd xb(21);
  Eigen::VectorXd u(6);
  Eigen::MatrixXd A, B;
  u << 9.81 / 0.16, 0, 0, 0.1, 0, 0;
  double h = 0.01;
  int N = 3.0 / h; // tf/h
  auto t0 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < N; ++i) {
    airm_system->Step(xb, 0, xa, u, h, 0, &A, &B);
    xa = xb;
  }
  double dt = std::chrono::duration<double>(
                  std::chrono::high_resolution_clock::now() - t0)
                  .count();
  std::cout << "Dt: " << dt << std::endl;
  // rp
  ASSERT_NEAR(xb(3), 0.2, 0.01);
  ASSERT_NEAR(xb(4), -0.2, 0.01);
  // rpydot
  ASSERT_NEAR(xb(9), 0.0, 0.01);
  ASSERT_NEAR(xb(10), 0.0, 0.01);
  ASSERT_NEAR(xb(11), 0.1, 0.01);
  // rpd
  ASSERT_DOUBLE_EQ(xb(12), 0.2);
  ASSERT_DOUBLE_EQ(xb(13), -0.2);
  // ja
  ASSERT_NEAR(xb(15), 0.5, 0.01);
  ASSERT_NEAR(xb(16), 0.5, 0.01);
  // jv
  ASSERT_NEAR(xb(17), 0.0, 0.01);
  ASSERT_NEAR(xb(18), 0.0, 0.01);
  // jad
  ASSERT_NEAR(xb(19), 0.5, 0.01);
  ASSERT_NEAR(xb(20), 0.5, 0.01);
  ASSERT_EQ(A.rows(), 21);
  ASSERT_EQ(A.cols(), 21);
  ASSERT_EQ(B.rows(), 21);
  ASSERT_EQ(B.cols(), 6);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
