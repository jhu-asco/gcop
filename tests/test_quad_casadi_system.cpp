#include "quad_casadi_system.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <memory>

#include <chrono>

using namespace gcop;

class TestQuadCasadiSystem : public testing::Test {
protected:
  virtual void SetUp() {
    parameters.resize(1);
    parameters << 0.16;
    kp_rpy<< 10, 10, 0;
    kd_rpy<< 5, 5, 2;
    quad_system.reset(new QuadCasadiSystem(parameters, kp_rpy, kd_rpy));
    quad_system->instantiateStepFunction();
  }
  template <class T> void assertVector(casadi::DM input, T expected_value) {
    std::vector<double> elems = input.get_elements();
    ASSERT_EQ(elems.size(), expected_value.size());
    ASSERT_GT(elems.size(), 0);
    for (int i = 0; i < elems.size(); ++i) {
      ASSERT_NEAR(elems[i], expected_value[i], 1e-7);
    }
  }
  std::unique_ptr<QuadCasadiSystem> quad_system;
  Eigen::VectorXd parameters;
  Eigen::Vector3d kp_rpy;
  Eigen::Vector3d kd_rpy;
};

TEST_F(TestQuadCasadiSystem, Initialize) {}

TEST_F(TestQuadCasadiSystem, InitializeCodeGen) {
  quad_system.reset(new QuadCasadiSystem(parameters, kp_rpy, kd_rpy, true));
  quad_system->instantiateStepFunction();
}

TEST_F(TestQuadCasadiSystem, TestBodyZAxes) {
  casadi::MX rpy_sym = casadi::MX::sym("rpy", 3);
  casadi::MX body_z_axis = quad_system->computeBodyZAxes(rpy_sym);
  casadi::Function z_axes_fun("body_z_axis_comp", {rpy_sym}, {body_z_axis});
  // flat
  casadi::DM rpy = std::vector<double>({0, 0, 1});
  casadi::DM z_axis = z_axes_fun(std::vector<casadi::DM>{rpy}).at(0);
  assertVector(z_axis, std::vector<double>{0, 0, 1});
  // Pitched
  rpy = std::vector<double>({0, M_PI / 2.0, 0});
  z_axis = z_axes_fun(std::vector<casadi::DM>{rpy}).at(0);
  assertVector(z_axis, std::vector<double>{1, 0, 0});
  // Rolled
  rpy = std::vector<double>({M_PI / 2.0, 0, 0});
  z_axis = z_axes_fun(std::vector<casadi::DM>{rpy}).at(0);
  assertVector(z_axis, std::vector<double>{0, -1, 0});
}

TEST_F(TestQuadCasadiSystem, TestStep) {
  // Use code generation
  quad_system.reset(new QuadCasadiSystem(parameters, kp_rpy, kd_rpy, true));
  quad_system->instantiateStepFunction();
  Eigen::VectorXd xa(15);
  // p,                rpy,         v,         rpydot,     rpyd
  xa << 0.1, 0.1, 0.1, 0, 0, 0, 0.1, -0.2, -0.1, 0, 0, 0, 0.2, -0.2, 0.5;
  Eigen::VectorXd xb(15);
  Eigen::VectorXd u(4);
  Eigen::MatrixXd A, B;
  u << 1, 0, 0, 0.1;
  double h = 0.01;
  int N = 3.0 / h; // tf/h
  auto t0 = std::chrono::high_resolution_clock::now();
  LoopTimer overall_loop_timer;
  for (int i = 0; i < N; ++i) {
    overall_loop_timer.loop_start();
    quad_system->Step(xb, 0, xa, u, h, 0, &A, &B);
    xa = xb;
    overall_loop_timer.loop_end();
  }
  std::cout << "Copy loop average time: "
            << (quad_system->copy_loop_timer_).average_loop_period()
            << std::endl;
  std::cout << "Function loop average time: "
            << (quad_system->fun_loop_timer_).average_loop_period()
            << std::endl;
  std::cout << "Overall loop average time: "
            << (overall_loop_timer).average_loop_period() << std::endl;
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
  ASSERT_EQ(A.rows(), 15);
  ASSERT_EQ(A.cols(), 15);
  ASSERT_EQ(B.rows(), 15);
  ASSERT_EQ(B.cols(), 4);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
