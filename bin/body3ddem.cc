#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "body3dview.h"
#include "utils.h"
#include "body3dcost.h"
#include "demview.h"
#include "pqpdem.h"
#include "constraintcost.h"
#include "multicost.h"
#include "unistd.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<Body3dState, 12, 6> Body3dDdp;
typedef PqpDem<Body3dState, 12, 6> Body3dPqpDem;
typedef ConstraintCost<Body3dState, 12, 6, Dynamic, 1> PqpDemCost;


void solver_process(Viewer* viewer)
{
  //  if (viewer)
  //    viewer->SetCamera(3.25, 48, -1.15, -1.35, -2.25);

  //  viewer->SetCamera(1, 23, -98.2, -77, -49);
  //  viewer->SetCamera(.87, 23, -57.2, -64, 3.7);
  //  Dem dem("maps/pic2.ppm", .25, 30);

  if (viewer)
    viewer->SetCamera(4.5, 45, -10, -8, -10);
  Dem dem("maps/simple.ppm", .5, 15);
  dem.Convolve(2);

  DemView dv(dem);
  if (viewer)
    viewer->Add(dv);

  int N = 64;      // discrete trajectory segments
  double tf = 15;   // time-horizon
  double h = tf/N; // time-step

  // system
  Body3d<> sys;

  // cost 
  Body3dState xf;
  xf.Clear();
  xf.p[0] = 5;  
  xf.p[1] = 5;  
  xf.p[2] = 3;

  Body3dCost<6> cost(sys, tf, xf);  

  double s = .01;
  cost.Qf(0,0) = 2; cost.Qf(1,1) = 2; cost.Qf(2,2) = 2;
  cost.Qf(3,3) = 20; cost.Qf(4,4) = 20; cost.Qf(5,5) = 20;

  cost.Qf(6,6) = 5; cost.Qf(7,7) = 5; cost.Qf(8,8) = 5;
  cost.Qf(9,9) = 50; cost.Qf(10,10) = 50; cost.Qf(11,11) = 50;

  cost.Qf = s*cost.Qf;
  
  cost.R(0,0) = .05; cost.R(1,1) = .05; cost.R(2,2) = .05; 
  cost.R(3,3) = .1; cost.R(4,4) = .1; cost.R(5,5) = .1;

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // states
  vector<Body3dState> xs(N+1);
  xs[0].Clear();
  xs[0].R.setIdentity();
  xs[0].p[0] = 20;
  xs[0].p[1] = 17;
  xs[0].p[2] = 10;

  // initial controls (e.g. hover at one place)
  vector<Vector6d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].setZero();
  }

  Body3dView<6> view(sys, &xs);
  if (viewer)
    viewer->Add(view);  

  double temp[50000];

  Body3dPqpDem pqp(dem, 1);

  double temp2[50000];

  PqpDemCost pqpcost(sys, tf, pqp);
  
  MultiCost<Body3dState, 12, 6> mcost(sys, tf);
  mcost.costs.push_back(&cost);   
  mcost.costs.push_back(&pqpcost);
    
  Body3dDdp ddp(sys, mcost, ts, xs, us);
  ddp.eps = 1e-3;
  ddp.mu = 10;

  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed

  getchar();

  for (int i = 0; i < 500; ++i) {
    timer_start(timer);
    ddp.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << ": took " << te << " us." << endl;        
    getchar();
    pqpcost.b = 2*pqpcost.b;
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
  viewer->frameName = "body3d/frames/frame";

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
