#include "aerial_manipulation_feedforward_system.h"
#include "load_eigen_matrix.h"
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
    ub.resize(6);
    ub << 1.2, 0.3, 0.3, 0.3, 0.7, 0.7;
    lb = -ub;
    airm_system.reset(new AerialManipulationFeedforwardSystem(
        parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, 0.7, lb, ub));
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
  Eigen::VectorXd lb;
  Eigen::VectorXd ub;
};

TEST_F(TestAerialManipulationFeedforwardSystem, Initialize) {}

TEST_F(TestAerialManipulationFeedforwardSystem, InitializeCodeGen) {
  airm_system.reset(new AerialManipulationFeedforwardSystem(
      parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, 0.7, lb, ub, true));
  airm_system->instantiateStepFunction();
}

TEST_F(TestAerialManipulationFeedforwardSystem, TestAgainstTensorflow) {
  // Get parameters states, controls, predicted values
  std::string folder_path =
      (std::string(DATA_PATH) + "/tensorflow_model_vars_16_8_tanh/");
  std::string model_out_path =
      (std::string(DATA_PATH) + "/model_output_no_residual_dynamics/");
  int index = 10;
  Eigen::VectorXd extended_state =
      loadEigenMatrix(model_out_path + "state_" + std::to_string(index));
  Eigen::VectorXd xa(21);
  for (int i = 0; i < 22; ++i) {
    if (i < 15) {
      xa[i] = extended_state[i];
    } else if (i > 15) { // Ignore i = 15 which is kt
      xa[i - 1] = extended_state[i];
    }
  }
  double kt = extended_state[15];
  parameters[0] = kt;
  Eigen::VectorXd kp_rp = loadEigenMatrix(folder_path + "/rpy_gains_kp_0");
  kp_rpy[0] = kp_rp[0];
  kp_rpy[1] = kp_rp[1];
  kp_rpy[2] = 0;
  kd_rpy = loadEigenMatrix(folder_path + "/rpy_gains_kd_0");
  kp_ja = loadEigenMatrix(folder_path + "/joint_gains_kp_0");
  kd_ja = loadEigenMatrix(folder_path + "/joint_gains_kd_0");
  // Constants
  double max_joint_vel = 0.7;
  airm_system.reset(new AerialManipulationFeedforwardSystem(
      parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, max_joint_vel, lb, ub, true));
  airm_system->instantiateStepFunction();
  // Load model data
  Eigen::MatrixXd model_predictions =
      loadEigenMatrix(model_out_path + "predictions_" + std::to_string(index));
  Eigen::MatrixXd controls =
      loadEigenMatrix(model_out_path + "controls_" + std::to_string(index));
  int N = model_predictions.rows();
  double h = 0.02;
  Eigen::Vector3d rpyd_prev(xa[12], xa[13], xa[14]);
  Eigen::Vector2d jad_prev(xa[19], xa[20]);

  Eigen::VectorXd xb(21);
  for (int i = 0; i < N; ++i) {
    Eigen::VectorXd u(6);
    Eigen::VectorXd u1 = controls.row(i);
    Eigen::Vector3d rpyd_curr(u1[0], u1[1], u1[2]);
    Eigen::Vector3d rpyd_dot = (rpyd_curr - rpyd_prev) / 0.02;
    rpyd_prev = rpyd_curr;
    Eigen::Vector2d jad_curr(u1[4], u1[5]);
    Eigen::Vector2d jad_dot = (jad_curr - jad_prev) / 0.02;
    jad_prev = jad_curr;
    u << (u1[3] * kt / 9.81), rpyd_dot[0], rpyd_dot[1], rpyd_dot[2], jad_dot[0],
        jad_dot[1];
    airm_system->Step(xb, 0, xa, u, h);
    Eigen::VectorXd predicted_sens = model_predictions.row(i);
    Eigen::VectorXd measured_sens(8);
    measured_sens << xb[0], xb[1], xb[2], xb[3], xb[4], xb[5], xb[15], xb[16];
    xa = xb;
    for (int i = 0; i < 8; ++i) {
      ASSERT_NEAR(predicted_sens[i], measured_sens[i], 1e-5);
    }
  }
}
TEST_F(TestAerialManipulationFeedforwardSystem, TestStep) {
  // Use code generation
  airm_system.reset(new AerialManipulationFeedforwardSystem(
      parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, 0.7, lb, ub, true));
  airm_system->instantiateStepFunction();
  Eigen::VectorXd xa(21);
  // p,                rpy,         v,         rpydot,     rpyd,          ja,       jv,        jad;
  xa << 0.1, 0.1, 0.1, 0, 0, 0, 0.1, -0.2, -0.1, 0, 0, 0, 0.2, -0.2, 0.5, 0.1, 0.2, -0.1, 0.0, 0.5, 0.5;
  Eigen::VectorXd xb(21);
  Eigen::VectorXd u(6);
  Eigen::MatrixXd A, B;
  u << 1, 0, 0, 0.1, 0, 0;
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
