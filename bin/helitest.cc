#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "heliview.h"
#include "utils.h"
#include "body3dcost.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<Body3dState, 12, 4> HeliDdp;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(3.25, 48, -1.15, -1.35, -2.25);

  int N = 64;      // discrete trajectory segments
  double tf = 5;   // time-horizon
  double h = tf/N; // time-step

  // system
  Heli sys;

  // cost 
  Body3dState xf;
  xf.Clear();

  Body3dCost<4> cost(sys, tf, xf);  

  cost.Qf(0,0) = 2; cost.Qf(1,1) = 2; cost.Qf(2,2) = 2;
  cost.Qf(3,3) = 50; cost.Qf(4,4) = 50; cost.Qf(5,5) = 50;

  cost.Qf(6,6) = 5; cost.Qf(7,7) = 5; cost.Qf(8,8) = 5;
  cost.Qf(9,9) = 150; cost.Qf(10,10) = 150; cost.Qf(11,11) = 150;

  cost.R(0,0) = .05; cost.R(1,1) = .05; cost.R(2,2) = .05; cost.R(3,3) = .1;

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // states
  vector<Body3dState> xs(N+1);
  xs[0].Clear();
  xs[0].p << -5, 0, 2;

  // initial controls (e.g. hover at one place)
  vector<Vector4d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].head(3).setZero();
    us[i][3] = 9.81*sys.m;
  }

  vector<Vector4d> usd = us;
  cost.SetReference(0, &usd);

  HeliDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1;

  HeliView view(sys, &ddp.xs);
  if (viewer)
    viewer->Add(view);  

  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed

  for (int i = 0; i < 20; ++i) {    
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
  viewer->frameName = "heli/frames/frame";

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
