#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "body2dview.h"
//#include "body2dcost.h"
#include "lqcost.h"
#include "utils.h"
#include "se2.h"
#include "params.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp< pair<Matrix3d, Vector3d>, 6, 2> Body2dDdp;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(-44.875, 35.875, 1.9, 2.450001, -11);

  int N = 32;
  double tf = 10;
  int iters = 50;

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);
  params.GetInt("iters", iters);

  double h = tf/N;

  SE2 &se2 = SE2::Instance();  
  Body2dForce<2> force;
  force.D << 1, .1, 5;
  params.GetVector3d("D", force.D);

  // B matrix
  double cr = .3;
  params.GetDouble("cr", cr);
  force.B << cr, -cr, 1, 1, 0, 0;

  Body2d<2> sys(&force);
  sys.d << 1.16, 0.72;

  params.GetVector3d("I", sys.I);

  params.GetVector2d("ulb", sys.U.lb);
  params.GetVector2d("uub", sys.U.ub);
  sys.U.bnd = true;

  Body2dState x0;
  VectorXd qv0(6);
  se2.q2g(x0.first, Vector3d(0, 2, 2));
  if (params.GetVectorXd("x0", qv0)) {
    se2.q2g(x0.first, qv0.head(3));
    x0.second = qv0.tail(3);
  }

  Body2dState xf;
  VectorXd qvf(6);
  if (params.GetVectorXd("xf", qvf)) {
    se2.q2g(xf.first, qvf.head(3));
    xf.second = qvf.tail(3);
  }

  // cost
  LqCost<Body2dState, 6, 2> cost(sys, tf, xf);

  //  Body2dCost cost(sys, tf, xf);  

  VectorXd Q(6);
  if (params.GetVectorXd("Q", Q))
    cost.Q = Q.asDiagonal();
  
  VectorXd Qf(6);
  if (params.GetVectorXd("Qf", Qf))
    cost.Qf = Qf.asDiagonal();
  
  VectorXd R(2);
  if (params.GetVectorXd("R", R)) 
    cost.R = R.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<pair<Matrix3d, Vector3d> > xs(N+1, x0);

  // controls
  Vector2d u(0,0);
  vector<Vector2d> us(N, u);
  for (int i = 0; i < N/2; ++i) {
    //    us[i] = Vector2d(0,1);
    //    us[N/2+i] = Vector2d(0,-1);
  }

  Body2dDdp ddp(sys, cost, ts, xs, us);

  
  ddp.mu = .01;
  params.GetDouble("mu", ddp.mu);
  params.GetDouble("eps", ddp.eps);

  Body2dView<2> view(sys, &ddp.xs);
  viewer->Add(view);

  struct timeval timer;

  ddp.debug = false; // turn off debug for speed

  getchar();

  for (int i = 0; i < iters; ++i) {

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
