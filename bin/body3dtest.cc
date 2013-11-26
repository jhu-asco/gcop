#include <iomanip>
#include <iostream>
#include "dmoc.h"
#include "viewer.h"
#include "body3dview.h"
#include "body3dcost.h"
#include "utils.h"
#include "so3.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<Body3dState, 12, 6> Body3dDmoc;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(3, 60, -2.55, -3.9, -5.05);

  int N = 32;
  double tf = 10;
  double h = tf/N;
  
  Body3d<> sys;

  Body3dState xf(Matrix3d::Identity(), Vector9d::Zero());

  //  Body3dCost<> cost(tf, xf);
  LqCost<Body3dState, 12, 6> cost(Body3dManifold::Instance(), tf, xf);

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Body3dState> xs(N+1);
  Vector3d e0(1.2, -2, 1);
  SO3::Instance().exp(xs[0].first, e0);
  xs[0].second[0] = 5;
  xs[0].second[1] = 5;
  xs[0].second[2] = 5;

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
  

  Body3dDmoc dmoc(sys, cost, ts, xs, us);
  dmoc.mu = 1;
  //  dmoc.a = .5;

  Body3dView<> view(sys, &dmoc.xs);
  viewer->Add(view);  

  struct timeval timer;

  //  dmoc.debug = false; // turn off debug for speed

  for (int i = 0; i < 20; ++i) {
    timer_start(timer);
    dmoc.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << " took: " << te << " us." << endl;    
    getchar();
  }

  //  for (int k = 0; k <= N; ++k)
  //    cout << dmoc.xs[k] << "|" << endl;  
  //  cout << "xf=" << dmoc.xs.back() << endl;
  
  
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
