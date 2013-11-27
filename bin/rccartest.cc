#include <iomanip>
#include <iostream>
#include "dmoc.h"
#include "viewer.h"
#include "rccarview.h"
#include "utils.h"
#include "rnlqcost.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<Vector4d, 4, 2> RccarDmoc;

void solver_process(Viewer* viewer)
{
  int N = 64;        // number of segments
  double tf = 20;    // time horizon
  double h = tf/N;   // time step

  Rccar sys;

  // final state
  Vector4d xf = Vector4d::Zero();

  // cost
  RnLqCost<4, 2> cost(tf, xf);
  cost.Q(0,0) = 0.1; cost.Q(1,1) = 0.1; cost.Q(2,2) = 0.1; cost.Q(3,3) = .5;
  cost.Qf(0,0) = 5; cost.Qf(1,1) = 5; cost.Qf(2,2) = 10; cost.Qf(3,3) = 1;
  cost.R(0,0) = 1; cost.R(1,1) = .1;  

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Vector4d> xs(N+1);
  // initial state
  xs[0][0] = .5; 
  xs[0][1] = .5;
  xs[0][2] = M_PI/4;
  xs[0][3] = 0;

  // initial controls
  vector<Vector2d> us(N);
  for (int i = 0; i < N/2; ++i) {
    us[i] = Vector2d(.01, .0);
    us[N/2+i] = Vector2d(-.01, .0);
  }
  
  RccarDmoc dmoc(sys, cost, ts, xs, us);  
  dmoc.mu = .1;

  RccarView view(sys, &dmoc.xs);
  
  viewer->Add(view);

  struct timeval timer;
  // dmoc.debug = false; // turn off debug for speed
  getchar();

  for (int i = 0; i < 30; ++i) {

    timer_start(timer);
    dmoc.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << " took: " << te << " us." << endl;
    getchar();
  }

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
