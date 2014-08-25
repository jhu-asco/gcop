#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "body2dview.h"
#include "body2dcost.h"
#include "utils.h"
#include "se2.h"
#include "params.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<pair<Matrix3d, Vector3d>, 6, 3> Body2dDdp;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(4.875, 33.875, 0.24999, -0.550001, -6);

  int N = 32;
  double tf = 10;
  double h = tf/N;
  int iters = 50;

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);
  params.GetInt("iters", iters);

  SE2 &se2 = SE2::Instance();  
  Body2dForce force;
  force.D(2)= 5;
  params.GetVector3d("D", force.D);

  Body2d sys(&force);

  M3V3d x0;
  VectorXd qv0(6);
  se2.q2g(x0.first, Vector3d(0, 2, 2));
  if (params.GetVectorXd("x0", qv0)) {
    se2.q2g(x0.first, qv0.head(3));
    x0.second = qv0.tail(3);
  }

  M3V3d xf;
  VectorXd qvf(6);
  if (params.GetVectorXd("xf", qvf)) {
    se2.q2g(xf.first, qvf.head(3));
    xf.second = qvf.tail(3);
  }

  // cost
  Body2dCost cost(sys, tf, xf);  

  VectorXd Q(6);
  if (params.GetVectorXd("Q", Q))
    cost.Q = Q.asDiagonal();
  
  VectorXd Qf(6);
  if (params.GetVectorXd("Qf", Qf))
    cost.Qf = Qf.asDiagonal();
  
  VectorXd R(3);
  if (params.GetVectorXd("R", R)) 
    cost.R = R.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<pair<Matrix3d, Vector3d> > xs(N+1, x0);

  // controls
  Vector3d u(0,0,0);
  vector<Vector3d> us(N, u);
  //  for (int i = 0; i < N/2; ++i) {
  //    us[i] = Vector3d(0.01,.1,0);
  //    us[N/2+i] = Vector3d(0.01,-.1,0);
  //  }  

  Body2dDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = .01;
  params.GetDouble("mu", ddp.mu);

  Body2dView view(sys, &ddp.xs);
  viewer->Add(view);

  struct timeval timer;

  ddp.debug = false; // turn off debug for speed

  getchar();

  for (int i = 0; i < 50; ++i) {

    timer_start(timer);
    ddp.Iterate();
    long te = timer_us(timer);

    cout << "Iteration #" << i << " took: " << te << " us." << endl;    
    getchar();  
  }

  //  for (int k = 0; k <= N; ++k)
  //    cout << ddp.xs[k] << "|" << endl;  
  //  cout << "xf=" << ddp.xs.back() << endl;
  
  
  cout << "done!" << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{

  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("../../bin/body2d.cfg");


#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "videos/sys";

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
