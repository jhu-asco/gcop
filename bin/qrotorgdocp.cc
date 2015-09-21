#include <iomanip>
#include <iostream>
#include "viewer.h"
#include "qrotorview.h"
#include "utils.h"
#include "body3dcost.h"
#include "params.h"
#include "glog/logging.h"


//#define USE_SGDOCP

#include "gdocp.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef GDocp<Body3dState, 12, 4> QrotorGDocp;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(43.25, 94, -1.15, -0.9, -3.75);


  int N = 32;      // discrete trajectory segments
  double tf = 2;   // time-horizon

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);  

  double h = tf/N; // time-step

  // system
  Qrotor sys;
  //  sys.symp = false;

  Body3dState x0;
  VectorXd qv0(12);
  params.GetVectorXd("x0", qv0);  
  SO3::Instance().q2g(x0.R, qv0.head(3));    
  x0.p = qv0.segment<3>(3); x0.w = qv0.segment<3>(6); x0.v = qv0.tail<3>(); 

  Body3dState xf;
  VectorXd qvf(12);
  params.GetVectorXd("xf", qvf);  
  SO3::Instance().q2g(xf.R, qvf.head(3));    
  xf.p = qvf.segment<3>(3); xf.w = qvf.segment<3>(6); xf.v = qvf.tail<3>(); 

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

  QrotorGDocp ocp(sys, cost, ts, xs, us);
  //  ocp.prec = true;
  ocp.mu = .01;

  QrotorView view(sys, &ocp.xs);
  if (viewer)
    viewer->Add(view);  

  struct timeval timer;
  //  gdocp.debug = false; // turn off debug for speed

  timer_start(timer);
  ocp.Iterate();
  long te = timer_us(timer);
  cout << "Solution" << ": took " << te << " us." << endl;        
  
  cout << "done!" << endl;
  cout << "Click on window and press 'a' to animate." << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{
  google::InitGoogleLogging(argv[0]);


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
