#include <iomanip>
#include <iostream>
#include "dmoc.h"
#include "viewer.h"
#include "hrotorview.h"
#include "utils.h"
#include "body3dcost.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<Body3dState, 12, 4> HrotorDmoc;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(3.25, 48, -1.15, -1.35, -2.25);

  int N = 64;      // discrete trajectory segments
  double tf = 5;   // time-horizon
  double h = tf/N; // time-step

  // system
  Hrotor sys;

  // cost 
  Body3dState xf(Matrix3d::Identity(), Vector9d::Zero());
  Body3dCost<4> cost(tf, xf);  
  cost.Qf(0,0) = 2; cost.Qf(1,1) = 2; cost.Qf(2,2) = 2;
  cost.Qf(3,3) = 20; cost.Qf(4,4) = 20; cost.Qf(5,5) = 20;

  cost.Qf(6,6) = 5; cost.Qf(7,7) = 5; cost.Qf(8,8) = 5;
  cost.Qf(9,9) = 50; cost.Qf(10,10) = 50; cost.Qf(11,11) = 50;

  cost.R(0,0) = .05; cost.R(1,1) = .05; cost.R(2,2) = .05; cost.R(3,3) = .1;

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // states
  vector<Body3dState> xs(N+1);
  xs[0].first.setIdentity();
  xs[0].second[0] = 5;
  xs[0].second[1] = 5;
  xs[0].second[2] = 2;

  // initial controls (e.g. hover at one place)
  vector<Vector4d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].head(3).setZero();
    us[i][3] = 9.81*sys.m;
  }
    
  HrotorDmoc dmoc(sys, cost, ts, xs, us);
  dmoc.mu = 1;

  HrotorView view(sys, &dmoc.xs);
  if (viewer)
    viewer->Add(view);  

  struct timeval timer;
  //  dmoc.debug = false; // turn off debug for speed

  for (int i = 0; i < 10; ++i) {    
    timer_start(timer);
    dmoc.Iterate();
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

#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "hrotor/frames/frame";

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
