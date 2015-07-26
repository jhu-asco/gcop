#include <iomanip>
#include <iostream>
#include "ddp.h"
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

typedef Ddp<Body3dState, 12, 4> HeliDdp;
typedef Ddp<Vector5d, 5, 2> CarDdp;


int N = 256;      // discrete trajectory segments
double tf = 8;   // time-horizon
double h = tf/N; // time-step


void car_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(-46, 71, 0.85, -0.75, -5.0);

  SE2 &se2 = SE2::Instance();  
  Car sys;

  Vector5d xf = Vector5d::Zero();
  RnLqCost<5, 2> cost(sys, tf, xf);
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

  CarDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 10;

  CarView view(sys, &ddp.xs);
  viewer->Add(view);

  struct timeval timer;

  ddp.debug = false; // turn off debug for speed

  for (int i = 0; i < 20; ++i) {
    ddp.Iterate();
  }

  //  for (int k = 0; k <= N; ++k)
  //    cout << ddp.xs[k] << "|" << endl;  
  //  cout << "xf=" << ddp.xs.back() << endl;
  while(1)
    usleep(10);    
}


void heli_process(Viewer* viewer)
{

  // system
  Heli sys;
  sys.symp = false;

  // cost 
  Body3dState xf;
  xf.Clear();
  xf.p << .5, 0, .4;

  Body3dCost<4> cost(sys, tf, xf);  
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
  xs[0].Clear();
  xs[0].p << -5, 0, 2;

  // initial controls (e.g. hover at one place)
  vector<Vector4d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].head(3).setZero();
    us[i][3] = 9.81*sys.m;
  }
    
  HeliDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1;

  HeliView view(sys, &ddp.xs);
  if (viewer)
    viewer->Add(view);  


  // replace if with while for RHC control
  //while(1) {    
  if(1) {    
    struct timeval timer;
    //  ddp.debug = false; // turn off debug for speed
    
    timer_start(timer);
    for (int i = 0; i < 20; ++i) {
      ddp.Iterate();
    }
    long te = timer_us(timer);
    cout << "Took " << te << " us." << endl;        
    // getchar();
    xs[0] = xs[1];

    std::rotate(us.begin(), us.begin() + 1, us.end());
    us.back() = us[us.size()-2];

    ddp.Update();
    cout << xs[0].R << " " << xs[0].p  << endl;
    cout << "Moved forward" << endl;
    //    getchar();
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
