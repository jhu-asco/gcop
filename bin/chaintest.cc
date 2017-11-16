#include <iomanip>
#include <iostream>
#include "viewer.h"
#include "chainview.h"
#include "utils.h"
#include "unistd.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Matrix<double, 8, 1> Vector8d;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(45.25, 37, -0.14, 0.05, -1.75);

  double h = .01;
  int N = 1000;     // discrete trajectory segments
  double tf = h*N;   // time-horizon


  // system
  Chain sys;
  sys.debug = true;

  MbsState x(4);
  x.gs[0].setIdentity();
  x.r[0] = 0;   
  x.r[1] = 0;
  x.r[2] = 0;
  sys.FK(x);

  // states
  vector<MbsState> xs(N+1, x);
  xs[0].vs[0].setZero();
  xs[0].dr[0] = .5;
  xs[0].dr[1] = .5;
  xs[0].dr[2] = .5;
  sys.KStep(xs[0], x, h);

  // initial controls (e.g. hover at one place)
  VectorXd u(8);
  u.setZero();
  vector<VectorXd> us(N, u);

  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed

  for (int i = 0; i < N; ++i) {    
    timer_start(timer);
    sys.Step(xs[i+1], i*h, xs[i], us[i], h);
    long te = timer_us(timer);
    cout << "Iteration #" << i << ": took " << te << " us." << endl;
  }
  cout << "dr=" << xs[1].dr << endl;

  cout << "done!" << endl;

  ChainView view(sys, &xs);
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
  viewer->frameName = "chain/frames/frame";

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
