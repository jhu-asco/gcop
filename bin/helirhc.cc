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
    viewer->SetCamera(6.75, 81, 2.7, -.9, -5.05);

  int N = 128;      // discrete trajectory segments
  double tf = 5;   // time-horizon
  double h = tf/N; // time-step

  // system
  Heli sys;
  sys.symp = false;

  // cost 
  Body3dState xf(Matrix3d::Identity(), Vector9d::Zero());
  Body3dCost<4> cost(sys, tf, xf);  
  cost.Qf(0,0) = .1; cost.Qf(1,1) = .1; cost.Qf(2,2) = .1;
  cost.Qf(3,3) = 5; cost.Qf(4,4) = 5; cost.Qf(5,5) = 5;

  cost.Qf(6,6) = 1; cost.Qf(7,7) = 1; cost.Qf(8,8) = 1;
  cost.Qf(9,9) = 10; cost.Qf(10,10) = 10; cost.Qf(11,11) = 10;

  cost.Q.setZero();

  cost.R(0,0) = 0.1; cost.R(1,1) = 0.1; cost.R(2,2) = 0.1; cost.R(3,3) = 1;


  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // states
  vector<Body3dState> xs(N+1);
  xs[0].first.setIdentity();
  xs[0].second[0] = -5;
  xs[0].second[1] = 0;
  xs[0].second[2] = 2;

  // initial controls (e.g. hover at one place)
  vector<Vector4d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].head(3).setZero();
    us[i][3] = 9.81*sys.m;
  }

  cost.SetReference(0, &us);
    
  HeliDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1;

  HeliView view(sys, &ddp.xs);
  if (viewer)
    viewer->Add(view);  


  while(1) {    
    struct timeval timer;
    //  ddp.debug = false; // turn off debug for speed
    
    timer_start(timer);
    for (int i = 0; i < 20; ++i) {
      ddp.Iterate();
    }
    long te = timer_us(timer);
    cout << "Took " << te << " us." << endl;        
    getchar();
    xs[0] = xs[1];
    //    xs[0].second[0] += .01;
    //    xs[0].second[1] += .01;
    xs[0].second[2] += .01;

    std::rotate(us.begin(), us.begin() + 1, us.end());
    us.back() = us[us.size()-2];

    ddp.Update();
    cout << xs[0].first << " " << xs[0].second  << endl;
    cout << "Moved forward" << endl;
    getchar();
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
