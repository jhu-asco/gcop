#include "airm_residual_network_model.h"
#include "load_eigen_matrix.h"
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
    double max_joint_vel = 0.7;
    int n_layers = 3;

    std::string folder_path =
        (std::string(DATA_PATH) + "/tensorflow_model_vars/");
    airm_model.reset(new AirmResidualNetworkModel(
        parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, max_joint_vel, n_layers,
        folder_path, Activation::tanh, true));
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

/*
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
}*/

TEST_F(TestAirmResidualNetworkModel, TestTrajectory) {
  // Get parameters states, controls, predicted values
  std::string folder_path =
      (std::string(DATA_PATH) + "/tensorflow_model_vars/");
  std::string model_out_path = (std::string(DATA_PATH) + "/model_output/");
  int index = 5;
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
  parameters[0] = extended_state[15];
  std::cout << "kt: " << parameters[0] << std::endl;
  Eigen::VectorXd kp_rp = loadEigenMatrix(folder_path + "/rpy_gains_kp_0");
  kp_rpy[0] = kp_rp[0];
  kp_rpy[1] = kp_rp[1];
  kp_rpy[2] = 0;
  kd_rpy = loadEigenMatrix(folder_path + "/rpy_gains_kd_0");
  kp_ja = loadEigenMatrix(folder_path + "/joint_gains_kp_0");
  kd_ja = loadEigenMatrix(folder_path + "/joint_gains_kd_0");
  // Constants
  double max_joint_vel = 0.7;
  int n_layers = 3;
  airm_model.reset(new AirmResidualNetworkModel(
      parameters, kp_rpy, kd_rpy, kp_ja, kd_ja, max_joint_vel, n_layers,
      folder_path, Activation::tanh, true));
  airm_model->instantiateStepFunction();
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
    u << u1[3], rpyd_dot[0], rpyd_dot[1], rpyd_dot[2], jad_dot[0], jad_dot[1];
    airm_model->Step(xb, 0, xa, u, h);
    // Verify xb is the same as prediction
    Eigen::VectorXd predicted_sens = model_predictions.row(i);
    Eigen::VectorXd measured_sens(8);
    measured_sens << xb[0], xb[1], xb[2], xb[3], xb[4], xb[5], xb[15], xb[16];
    /*
    std::cout<<"i:"<<i<<std::endl;
    std::cout<<"predicted: "<<predicted_sens.transpose()<<std::endl;
    std::cout<<"observed: "<<measured_sens.transpose()<<std::endl;
    //std::cout<<"observed: "<<xb.transpose()<<std::endl;
    std::cout<<"u1: "<<u1.transpose()<<std::endl;
    for(int i = 0; i < 8; ++i) {
         ASSERT_NEAR(predicted_sens[i], measured_sens[i], 1e-2);
     }*/
    xa = xb;
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
