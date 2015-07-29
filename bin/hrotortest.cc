#include <iomanip>
#include <iostream>
#include "viewer.h"
#include "hrotorview.h"
#include "utils.h"
#include "body3dcost.h"

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
typedef SDdp<Body3dState, 12, 4> HrotorDdp;
#else
typedef Ddp<Body3dState, 12, 4> HrotorDdp;
#endif

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(2.75, 60.75, -1.7, -1.4, -5.3);

  int N = 64;      // discrete trajectory segments
  double tf = 5;   // time-horizon
  double h = tf/N; // time-step

  // system
  Hrotor sys;

  // cost 
  Body3dState xf;
  xf.Clear();

  Body3dCost<4> cost(sys, tf, xf);  
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
  xs[0].Clear();
  xs[0].p << 5, 5, 2;

  // initial controls (e.g. hover at one place)
  vector<Vector4d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].head(3).setZero();
    us[i][3] = 9.81*sys.m;
  }
    
  HrotorDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1;

  HrotorView view(sys, &ddp.xs);
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
