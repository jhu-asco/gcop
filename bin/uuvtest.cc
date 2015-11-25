#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "uuvview.h"
#include "utils.h"
#include "se3.h"
#include "lqcost.h"
#include "params.h"

#define DISP

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<UuvState, 12, 6> UuvDdp;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(42, 65, -1.55, -2.45, -9);

  int N = 64;        // number of segments
  double tf = 10;    // time horizon

  int iters = 30;    // DDP iterations

  Uuv<> sys;         // default UUV system

  // see config file, e.g. jhurov.cfg for description of the parameters
  params.GetInt("N", N);  
  params.GetDouble("tf", tf);
  params.GetInt("iters", iters);

  params.GetMatrix6d("M", sys.I);
  params.GetVector6d("H", sys.d);
  params.GetVector3d("b", sys.b);
  params.GetVector3d("g", sys.fp);

  params.GetVector6d("ulb", sys.U.lb);
  params.GetVector6d("uub", sys.U.ub);
  sys.U.bnd = true;

  // create underactuation
  //  sys.U.lb[4] = -.0001;
  //  sys.U.ub[4] = .0001;

  // time-step
  double h = tf/N;

  // initial state
  UuvState x0;
  VectorXd qv0(12);
  if (params.GetVectorXd("x0", qv0)) {
    SE3::Instance().q2g(x0.g, qv0.head<6>());
    x0.v = qv0.tail<6>();
  }

  // final state
  UuvState xf;
  VectorXd qvf(12);
  if (params.GetVectorXd("xf", qvf)) {
    SE3::Instance().q2g(xf.g, qvf.head<6>());
    xf.v = qvf.tail<6>();
  }


  //  linear quadratic cost
  LqCost<UuvState, 12, 6> cost(sys, tf, xf);

  VectorXd Q(12);
  if (params.GetVectorXd("Q", Q))
    cost.Q = Q.asDiagonal();
  
  VectorXd Qf(12);
  if (params.GetVectorXd("Qf", Qf))
    cost.Qf = Qf.asDiagonal();
  
  VectorXd R(6);
  if (params.GetVectorXd("R", R)) 
    cost.R = R.asDiagonal();


  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<UuvState> xs(N+1);
  xs[0] = x0;

  // initial controls, e.g. set to zero
  vector<Vector6d> us(N);
  for (int i = 0; i < N/2; ++i) {
    for (int j = 0; j < 3; ++j) {
      us[i][j] = .0;
    }
    for (int j = 3; j < 6; ++j) {
      us[i][j] = 0;
    }
    us[N/2+i] = -us[i];
  }
  

  UuvDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = .001;
  //  ddp.a = .5;
  
#ifdef DISP
  UuvView<> view(sys, &ddp.xs);
  viewer->Add(view);  
#endif
  struct timeval timer;
  //getchar();

  //  ddp.debug = false; // turn off debug for speed

  for (int i = 0; i < iters; ++i) {
    timer_start(timer);
    ddp.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << " took: " << te << " us." << endl;    
  }

  //  for (int k = 0; k <= N; ++k)
  //    cout << ddp.xs[k] << "|" << endl;  
  //  cout << "xf=" << ddp.xs.back() << endl;
  
  
  cout << "done!" << endl;
  while(1)
    usleep(10);    
}




int main(int argc, char** argv)
{

  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("../../bin/jhurov.cfg"); 


#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "../logs/uuv/frames/frame";
  viewer->displayName = "../logs/uuv/display/display";

  pthread_t dummy;
  pthread_create( &dummy, NULL, (void *(*) (void *)) solver_process, viewer);

#else
  solver_process(0);
#endif


#ifdef DISP
  viewer->Start();
#endif


  return 0;
}
