#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "gcarview.h"
#include "utils.h"
#include "se2.h"
#include "lqcost.h"
#include "unistd.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<GcarState, 4, 2> GcarDdp;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(4.875, 33.875, 0.24999, -0.550001, -6);

  int N = 32;
  double tf = 10;
  double h = tf/N;

  SE2 &se2 = SE2::Instance();  
  Gcar sys;

  GcarState xf;

  LqCost< GcarState, 4, 2> cost(sys, tf, xf);
  cost.Q.setZero();
  cost.R(0,0) = .005;
  cost.R(1,1) = .001;

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<GcarState> xs(N+1);
  Vector3d q0(0, 2, -1);
  se2.q2g(xs[0].g, q0); // set initial pose

  // controls
  vector<Vector2d> us(N);
  for (int i = 0; i < N/2; ++i) {
    us[i].setZero();
    //    us[i] = Vector2d(0.01,.1);
    //    us[N/2+i] = Vector2d(0.01,-.1);
  }
  
  GcarDdp ddp(sys, cost, ts, xs, us);
  //  ddp.mu = 1;

  GcarView view(sys, &ddp.xs, &ddp.us);
  viewer->Add(view);

  struct timeval timer;

  ddp.debug = false; // turn off debug for speed

  getchar();

  for (int i = 0; i < 50; ++i) {

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
