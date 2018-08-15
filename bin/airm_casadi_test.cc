#include "aerial_manipulation_feedforward_system.h"
#include "airm_residual_network_model.h"
#include "ddp.h"
#include "load_eigen_matrix.h"
#include "loop_timer.h"
#include "lqcost.h"
#include <Eigen/Dense>
#include <iostream>

#define USE_NN_MODEL

using namespace std;
using namespace Eigen;
using namespace gcop;

#ifdef USE_NN_MODEL
using System = AirmResidualNetworkModel;
#else
using System = AerialManipulationFeedforwardSystem;
#endif

typedef Ddp<VectorXd> AirmDdp;

void solver_process() {

  int N = 100;   // number of segments
  double tf = 2; // time horizon

  int iters = 50;
  bool use_code_generation = true;

  double h = tf / N;                     // time step
  Eigen::VectorXd dynamic_parameters(1); // kt
  dynamic_parameters << 0.15;
  Eigen::Vector3d kp_rpy;
  kp_rpy << 10, 10, 0;
  Eigen::Vector3d kd_rpy;
  kd_rpy << 5, 5, 2;
  Eigen::Vector2d kp_ja;
  kp_ja << 10, 10;
  Eigen::Vector2d kd_ja;
  kd_ja << 5, 5;
  Eigen::VectorXd ub(6);
  ub << 1.2, 0.6, 0.6, 0.6, 0.7, 0.7;
  Eigen::VectorXd lb = -ub;
  lb[0] = 0.8;
#ifdef USE_NN_MODEL
  std::string folder_path =
      (std::string(DATA_PATH) + "/tensorflow_model_vars_16_8_tanh/");
  Eigen::VectorXd kp_rp = loadEigenMatrix(folder_path + "/rpy_gains_kp_0");
  kp_rpy[0] = kp_rp[0];
  kp_rpy[1] = kp_rp[1];
  kp_rpy[2] = 0;
  kd_rpy = loadEigenMatrix(folder_path + "/rpy_gains_kd_0");
  kp_ja = loadEigenMatrix(folder_path + "/joint_gains_kp_0");
  kd_ja = loadEigenMatrix(folder_path + "/joint_gains_kd_0");
  double max_joint_vel = 0.7;
  int n_layers = 3;
  AirmResidualNetworkModel sys(dynamic_parameters, kp_rpy, kd_rpy, kp_ja, kd_ja,
                               max_joint_vel, n_layers, folder_path, lb, ub,
                               Activation::tanh, use_code_generation);
#else
  AerialManipulationFeedforwardSystem sys(dynamic_parameters, kp_rpy, kd_rpy,
                                          kp_ja, kd_ja, 0.7, lb, ub,
                                          use_code_generation);
#endif

  sys.instantiateStepFunction();

  // initial state
  VectorXd x0(21);
  // xyz,  rpy,  vxyz,  rpydot,  rpyd, ja, jv, jad
  x0 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1;

  // final state
  VectorXd xf(21);
  // xyz,  rpy,  vxyz,  rpydot,  rpyd, ja, jv, jad
  xf << 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1;

  // cost
  LqCost<VectorXd> cost(sys, tf, xf);
  VectorXd Q(21);
  Q << 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0;
  Q = 0.1 * Q;
  cost.Q = Q.asDiagonal();

  VectorXd Qf(21);
  Qf << 100, 100, 100, 800, 800, 800, 100, 100, 100, 100, 100, 100, 0.1, 0.1,
      0.1, 400, 400, 100, 100, 400, 400;

  cost.Qf = Qf.asDiagonal();

  VectorXd R(6);
  R << 6, 4.0, 4.0, 4.0, 4.0, 4.0;
  cost.R = R.asDiagonal();

  LoopTimer loop_timer;

  // times
  vector<double> ts(N + 1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k * h;

  // states
  vector<VectorXd> xs(N + 1);
  // initial state
  xs[0] = x0;

  // initial controls
  vector<VectorXd> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].resize(6);
    us[i] << 1.0, 0, 0, 0, 0, 0;
  }

  AirmDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1e-3;

  ddp.debug = false; // turn off debug for speed

  for (int i = 0; i < iters; ++i) {
    loop_timer.loop_start();
    ddp.Iterate();
    loop_timer.loop_end();
    cout << "Iters[" << (i + 1) << "]: Cost: [" << ddp.J << "]: dT: ["
         << loop_timer.average_loop_period() << "]" << endl;
  }
  for (int i = 0; i < N; ++i) {
    cout << "Index: " << i << "U: " << us[i].transpose() << std::endl;
  }

  cout << "Final state: " << endl;
  cout << xs[N].transpose() << endl;

  cout << "done!" << endl;
}

#define DISP

int main(int argc, char **argv) {

  solver_process();

  return 0;
}
