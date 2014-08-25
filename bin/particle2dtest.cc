#include <iostream>
#include "ddp.h"
#include "particle2d.h"
#include "rnlqcost.h"
#include "viewer.h"
#include "particle2dview.h"
#include "utils.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<Vector4d, 4, 2> Particle2dDdp;

void solver_process(Viewer* viewer)
{
  int N = 32;
  double tf = 10;
  double h = tf/N;

  Particle2d sys;

  Vector4d xf = Vector4d::Zero();
  RnLqCost<4, 2> cost(sys, tf, xf);
  cost.Q = Matrix4d(Vector4d(.01, .01, .005, .005).asDiagonal());
  cost.Qf = Matrix4d(Vector4d(1, 1, 5, 5).asDiagonal());
  cost.R = Matrix2d(Vector2d(.1, .1).asDiagonal());  


  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Vector4d> xs(N+1);
  xs[0] << -1, -1, .1, 0;

  // initial guess for the controls (almost any guess should work for this system)
  vector<Vector2d> us(N);
  for (int i = 0; i < N/2; ++i) {
    us[i] = Vector2d(.05, .025);  
    us[N/2+i] = Vector2d(-.05, -.025);
  }
  

  Particle2dDdp ddp(sys, cost, ts, xs, us);

  // should converge in one iteration since it is a linear system
  int iters = 5;  

  Particle2dView view(sys, &ddp.xs);
  viewer->Add(view);  

  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed

  for (int i = 0; i < iters; ++i) {

    timer_start(timer);
    ddp.Iterate();
    long te = timer_us(timer);

    cout << "Iteration #" << i << " took: " << te << " us." << endl;
    getchar();
  }
  //  for (int k = 0; k <= N; ++k)
  //cout << "xf=" << endl;
  //  cout << ddp.xs.back() << endl;

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
