#include "lqcost.h"
#include "params.h"
#include "quad_casadi_system.h"
#include "unistd.h"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

//#define USE_SDDP

#ifdef USE_SDDP
#include "sddp.h"
#else
#include "ddp.h"
#endif

using namespace std;
using namespace Eigen;
using namespace gcop;

#ifdef USE_SDDP
typedef SDdp<VectorXd> QRotorDdp;
#else
typedef Ddp<VectorXd> QRotorDdp;
#endif

// Params params;

void solver_process() {

  int N = 100;   // number of segments
  double tf = 1; // time horizon

  int iters = 50;
  bool use_code_generation = true;

  /*params.GetInt("N", N);
  params.GetDouble("tf", tf);

  params.GetInt("iters", iters);
  */

  double h = tf / N; // time step
  Eigen::VectorXd sys_params(7); // kt, kp, kd
  sys_params << 0.16, 10, 10, 10, 5, 5, 2;

  QuadCasadiSystem sys(sys_params, use_code_generation);
  sys.instantiateStepFunction();

  //  sys.U.lb[1] = tan(-M_PI/5);
  //  sys.U.ub[1] = tan(M_PI/5);

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
  Q.setZero();
  // Q.segment<3>(3) = Vector3d::Ones();
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
    us[i] << 9.81 / 0.16, 0, 0, 0;
  }

  QRotorDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1e-3;
#ifdef USE_SDDP
  ddp.duscale.resize(4);
  ddp.duscale << 0.01, 0.1, 0.1, 0.1;
#endif
  // ddp.mu = 1e-3;

  // struct timeval timer;
  // ddp.debug = false; // turn off debug for speed
  //  getchar();

  // timer_start(timer);
  for (int i = 0; i < iters; ++i) {
    ddp.Iterate();
    //    getchar();
  }

  // long te = timer_us(timer);
  // cout << "Iterations" << iters << " took: " << te << " us." << endl;
  //  getchar();

  cout << "Final state: " << endl;
  cout << xs[N].transpose() << endl;

  //  xs[1][3]  velocity
  // atan(us[0][1]) steering angle

  cout << "done!" << endl;
}

#define DISP

int main(int argc, char **argv) {

  /*if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("../../bin/rccar.cfg");
    */

  solver_process();

  return 0;
}
