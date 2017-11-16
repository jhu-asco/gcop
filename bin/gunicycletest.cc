#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "gunicycleview.h"
#include "gunicyclecost.h"
#include "utils.h"
#include "se2.h"
#include "unistd.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<pair<Matrix3d, Vector2d>, 5, 2> GunicycleDdp;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(4.875, 33.875, 0.24999, -0.550001, -6);

  int N = 32;
  double tf = 10;
  double h = tf/N;

  SE2 &se2 = SE2::Instance();  
  Gunicycle sys;

  M3V2d xf(se2.Id, Vector2d::Zero());
  GunicycleCost cost(sys, tf, xf);  

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<pair<Matrix3d, Vector2d> > xs(N+1);
  Vector3d q0(0, 0, -2);
  se2.q2g(xs[0].first, q0);
  xs[0].second[0] = 0;
  xs[0].second[1] = 0;

  // controls
  vector<Vector2d> us(N);
  for (int i = 0; i < N/2; ++i) {
    us[i] = Vector2d(0.01,.1);
    us[N/2+i] = Vector2d(0.01,-.1);
  }
  

  GunicycleDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 10;

  GunicycleView view(sys, &ddp.xs);
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
