#include <iomanip>
#include <iostream>
#include "viewer.h"
#include "body3dview.h"
#include "body3dcost.h"
#include "utils.h"
#include "so3.h"
#include "unistd.h"

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
typedef SDdp<Body3dState, 12, 6> Body3dDdp;
#else
typedef Ddp<Body3dState, 12, 6> Body3dDdp;
#endif

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(-2.5, 81, -1.8, -2.15, -6.3);

  int N = 32;
  double tf = 10;
  double h = tf/N;
  
  Body3d<> sys;

  Body3dState xf;
  xf.Clear();

  //  Body3dCost<> cost(tf, xf);
  LqCost<Body3dState, 12, 6> cost(sys, tf, xf);

  // times

  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Body3dState> xs(N+1);
  xs[0].Clear();
  Vector3d e0(1.2, -2, 1);
  SO3::Instance().exp(xs[0].R, e0);
  xs[0].p << 5,5,5;

  // controls 
  vector<Vector6d> us(N);
  for (int i = 0; i < N/2; ++i) {
    for (int j = 0; j < 3; ++j) {
      us[i][j] = .0;
    }
    for (int j = 3; j < 6; ++j) {
      us[i][j] = 0;
    }
    us[N/2+i] = -us[i];
  }
  

  Body3dDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1;
  //  ddp.a = .5;

  Body3dView<> view(sys, &ddp.xs);
  viewer->Add(view);  

  struct timeval timer;

  //  ddp.debug = false; // turn off debug for speed

  for (int i = 0; i < 20; ++i) {
    timer_start(timer);
    ddp.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << " took: " << te << " us." << endl;    
    //    getchar();
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
