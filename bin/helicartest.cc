#include <iomanip>
#include <iostream>
#include "dmoc.h"
#include "viewer.h"
#include "heliview.h"
#include "utils.h"
#include "body3dcost.h"
#include "carview.h"
#include "utils.h"
#include "se2.h"
#include "rnlqcost.h"


using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<Body3dState, 12, 4> HeliDmoc;
typedef Dmoc<Vector5d, 5, 2> CarDmoc;


int N = 256;      // discrete trajectory segments
double tf = 8;   // time-horizon
double h = tf/N; // time-step


void car_process(Viewer* viewer)
{
  //  if (viewer)
  // viewer->SetCamera(4.875, 33.875, 0.24999, -0.550001, -6);

  SE2 &se2 = SE2::Instance();  
  Car sys;

  Vector5d xf = Vector5d::Zero();
  RnLqCost<5, 2> cost(tf, xf);
  cost.Q(0,0) = 0.1; cost.Q(1,1) = 0.1; cost.Q(2,2) = 0.1; cost.Q(3,3) = .5; cost.Q(4,4) = .5;
  cost.Qf(0,0) = 10; cost.Qf(1,1) = 10; cost.Qf(2,2) = 1; cost.Qf(3,3) = 1; cost.Qf(4,4) = 0;  
  cost.R(0,0) = .01; cost.R(1,1) = .001;  

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Vector5d>  xs(N+1);
  xs[0][0] = -10; 
  xs[0][1] = 5;
  xs[0][2] = 0;
  xs[0][3] = 0;
  xs[0][4] = 0;

  // controls
  vector<Vector2d> us(N);
  for (int i = 0; i < N/2; ++i) {
    us[i] = Vector2d(0,0);
    us[N/2+i] = Vector2d(0,0);
  }

  CarDmoc dmoc(sys, cost, ts, xs, us);
  dmoc.mu = 10;

  CarView view(sys, &dmoc.xs);
  viewer->Add(view);

  struct timeval timer;

  dmoc.debug = false; // turn off debug for speed

  for (int i = 0; i < 20; ++i) {
    dmoc.Iterate();
  }

  //  for (int k = 0; k <= N; ++k)
  //    cout << dmoc.xs[k] << "|" << endl;  
  //  cout << "xf=" << dmoc.xs.back() << endl;
  while(1)
    usleep(10);    
}


void heli_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(6.75, 81, 2.7, -.9, -5.05);

  // system
  Heli sys;
  sys.symp = false;

  // cost 
  Body3dState xf(Matrix3d::Identity(), Vector9d::Zero());
  xf.second[0] = .5;
  xf.second[2] = .4;
  Body3dCost<4> cost(tf, xf);  
  cost.Qf(0,0) = .1; cost.Qf(1,1) = .1; cost.Qf(2,2) = .1;
  cost.Qf(3,3) = 5; cost.Qf(4,4) = 5; cost.Qf(5,5) = 5;

  cost.Qf(6,6) = 1; cost.Qf(7,7) = 1; cost.Qf(8,8) = 1;
  cost.Qf(9,9) = 10; cost.Qf(10,10) = 10; cost.Qf(11,11) = 10;

  cost.R(0,0) = 0; cost.R(1,1) = 0; cost.R(2,2) = 0; cost.R(3,3) = 0;

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
    
  HeliDmoc dmoc(sys, cost, ts, xs, us);
  dmoc.mu = 1;

  HeliView view(sys, &dmoc.xs);
  if (viewer)
    viewer->Add(view);  


  while(1) {    
    struct timeval timer;
    //  dmoc.debug = false; // turn off debug for speed
    
    timer_start(timer);
    for (int i = 0; i < 20; ++i) {
      dmoc.Iterate();
    }
    long te = timer_us(timer);
    cout << "Took " << te << " us." << endl;        
    getchar();
    xs[0] = xs[1];

    std::rotate(us.begin(), us.begin() + 1, us.end());
    us.back() = us[us.size()-2];

    dmoc.Update();
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
  pthread_create( &dummy, NULL, (void *(*) (void *)) car_process, viewer);


  pthread_t dummy2;
  pthread_create( &dummy2, NULL, (void *(*) (void *)) heli_process, viewer);

#else
  pthread_t dummy;
  pthread_create( &dummy, NULL, (void *(*) (void *)) car_process, 0);


  pthread_t dummy2;
  pthread_create( &dummy2, NULL, (void *(*) (void *)) heli_process, 0);

#endif


#ifdef DISP
  viewer->Start();
#endif


  return 0;
}
