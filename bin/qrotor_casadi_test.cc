#include "ddp.h"
#include "lqcost.h"
#include "quad_casadi_system.h"
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<VectorXd> QRotorDdp;

void solver_process() {

  int N = 100;   // number of segments
  double tf = 1; // time horizon

  int iters = 50;
  bool use_code_generation = true;

  double h = tf / N; // time step
  Eigen::VectorXd sys_params(1); // kt, kp, kd
  sys_params << 0.16;
  Eigen::Vector3d kp_rpy;
  kp_rpy<<10, 10, 10;
  Eigen::Vector3d kd_rpy;
  kd_rpy<<5, 5, 2;

  QuadCasadiSystem sys(sys_params, kp_rpy, kd_rpy, use_code_generation);
  sys.instantiateStepFunction();

  // initial state
  VectorXd x0(15);
  // xyz,  rpy,  vxyz,  rpydot,  rpyd
  x0 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  // final state
  VectorXd xf(15);
  // xyz,  rpy,  vxyz,  rpydot,  rpyd
  xf << 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  // cost
  LqCost<VectorXd> cost(sys, tf, xf);
  VectorXd Q(15);
  Q << 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
  cost.Q = Q.asDiagonal();

  VectorXd Qf(15);
  Qf << 100, 100, 100, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1;

  cost.Qf = Qf.asDiagonal();

  VectorXd R(4);
  R << 1e-4, 0.1, 0.1, 0.1;
  cost.R = R.asDiagonal();

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
    us[i].resize(4);
    us[i] << 1.0, 0, 0, 0;
  }

  QRotorDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1e-3;

  ddp.debug = false; // turn off debug for speed

  for (int i = 0; i < iters; ++i) {
    ddp.Iterate();
    cout << "Iters[" << (i + 1) << "]: Cost: " << ddp.J << endl;
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
