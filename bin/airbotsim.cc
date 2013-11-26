#include <iomanip>
#include <iostream>
#include "viewer.h"
#include "airbotview.h"
#include "utils.h"
#include "cylinderview.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Matrix<double, 8, 1> Vector8d;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(45.25, 37, -0.14, 0.05, -1.75);

  double h = .005;
  int N = 2;     // discrete trajectory segments
  double tf = h*N;   // time-horizon

  // system
  Airbot sys;

  MbsState x(7);
  x.gs[0].setIdentity();
  x.r[0] = 0;   
  x.r[1] = 0;
  sys.FK(x);

  // states
  vector<MbsState> xs(N+1, x);
  xs[0].vs[0].setZero();
  xs[0].dr[0] = 0;
  xs[0].dr[1] = 0;
  sys.Rec(xs[0], h);

  //  CylinderView cylv;
  //  SE3::Instance().rpyxyz2g(cylv.g0, Vector3d(0, 0, 0), Vector3d(0, 0, -sys.links[i+1].ds[0]/2));
  //  viewer->Add(cylv);

  // initial controls (e.g. hover at one place)
  VectorXd u(10);
  u.setZero();
  vector<VectorXd> us(N, u);
  //  us[0][4] = .1;
  //  us[1][4] = .1;
  //  us[2][4] = .1;

  struct timeval timer;
  //  dmoc.debug = false; // turn off debug for speed

  for (int i = 0; i < N; ++i) {    
    timer_start(timer);
    sys.Step(xs[i+1], i*h, xs[i], us[i], h);
    long te = timer_us(timer);
    cout << "Iteration #" << i << ": took " << te << " us." << endl;
  }
  cout << "dr=" << xs[1].dr << endl;

  cout << "done!" << endl;

  AirbotView view(sys, &xs);
  if (viewer)
    viewer->Add(view);  


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
  viewer->frameName = "../../logs/airbot/frames/frame";

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
