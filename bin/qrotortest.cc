#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "qrotorview.h"
#include "utils.h"
#include "body3dcost.h"
#include "params.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<Body3dState, 12, 4> QrotorDdp;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(3.25, 48, -1.15, -1.35, -2.25);

  int N = 32;      // discrete trajectory segments
  double tf = 2;   // time-horizon

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);  

  double h = tf/N; // time-step

  // system
  Qrotor sys;

  // cost 
  VectorXd qv0(12);
  params.GetVectorXd("x0", qv0);

  Body3dState x0;
  SO3::Instance().q2g(x0.first, qv0.head(3));
  x0.second= qv0.tail(9);

  VectorXd qvf(12);
  params.GetVectorXd("xf", qvf);  

  Body3dState xf;
  SO3::Instance().q2g(xf.first, qvf.head(3));
  xf.second= qvf.tail(9);

  Body3dCost<4> cost(sys, tf, xf);  

  VectorXd Q(12);
  VectorXd R(4);  
  VectorXd Qf(12);

  params.GetVectorXd("Q", Q);  
  params.GetVectorXd("R", R);
  params.GetVectorXd("Qf", Qf);

  cost.Q = Q.asDiagonal();
  cost.R = R.asDiagonal();
  cost.Qf = Qf.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // states
  vector<Body3dState> xs(N+1, x0);

  // initial controls (e.g. hover at one place)
  vector<Vector4d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].head(3).setZero();
    us[i][3] = 9.81*sys.m;
  }
    
  QrotorDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1;

  QrotorView view(sys, &ddp.xs);
  if (viewer)
    viewer->Add(view);  

  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed

  for (int i = 0; i < 10; ++i) {    
    timer_start(timer);
    ddp.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << ": took " << te << " us." << endl;        
  }

  cout << "done!" << endl;
  cout << "Click on window and press 'a' to animate." << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{
  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("qrotor.cfg");


#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "qrotor/frames/frame";

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
