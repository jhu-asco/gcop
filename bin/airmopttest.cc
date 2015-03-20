#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "airmview.h"
#include "utils.h"
#include "lqcost.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<MbsState> AirmDdp;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(-54, 63, 0.26, 0.2, -1.4);

  int N = 128;      // discrete trajectory segments
  double tf = 5;   // time-horizon
  double h = tf/N; // time-step


  // system
  Airm sys;

  MbsState xf(3);
  xf.gs[0].setIdentity();
  xf.r[0] = -1;
  xf.r[1] = -1;
  xf.vs[0].setZero();
  xf.dr.setZero();

  LqCost<MbsState> cost(sys, tf, xf);

  //  cost.Qf(3,3) = 50; cost.Qf(4,4) = 50; cost.Qf(5,5) = 50;
  //  cost.Qf(11,11) = 5; cost.Qf(12,12) = 5; cost.Qf(13,13) = 5;
  
  cost.Qf(0,0) = 2; cost.Qf(1,1) = 2; cost.Qf(2,2) = 2;
  cost.Qf(3,3) = 20; cost.Qf(4,4) = 20; cost.Qf(5,5) = 20;
  cost.Qf(6,6) = 10; cost.Qf(7,7) = 10;

  cost.Qf(8,8) = 5; cost.Qf(9,9) = 5; cost.Qf(10,10) = 5;
  cost.Qf(11,11) = 20; cost.Qf(12,12) = 20; cost.Qf(13,13) = 20;
  cost.Qf(14,14) = 50; cost.Qf(15,15) = 50;   
  

  cost.Q(0,0) = 2; cost.Q(1,1) = 2; cost.Q(2,2) = 2;
  cost.Q(3,3) = 20; cost.Q(4,4) = 20; cost.Q(5,5) = 20;
  cost.Q(6,6) = 1; cost.Q(7,7) = 1; 

  cost.Q(8,8) = 5; cost.Q(9,9) = 5; cost.Q(10,10) = 5;
  cost.Q(11,11) = 20; cost.Q(12,12) = 20; cost.Q(13,13) = 20;
  cost.Q(14,14) = 10; cost.Q(15,15) = 10;

  cost.R(0,0) = .01; cost.R(1,1) = .01; cost.R(2,2) = .01; cost.R(3,3) = .0;
  cost.R(4,4) = .0;  cost.R(5,5) = .0;

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  MbsState x(3);
  SE3::Instance().rpyxyz2g(x.gs[0], Vector3d(0,0,0), Vector3d(-2,0,0));
  x.r[0] = 0;   
  x.r[1] = 0;
  x.vs[0].setZero();
  x.dr.setZero();  
  sys.Rec(x,h);

  // states
  vector<MbsState> xs(N+1, x);
  //  xs[0].vs[0].setZero();
  //  xs[0].dr[0] = 0;
  //  xs[0].dr[1] = 0;
  //  sys.KStep(xs[0], x, h);

  double m = sys.links[0].m + sys.links[1].m  + sys.links[2].m;
  
  // initial controls (e.g. hover at one place)
  VectorXd u(6);
  u.setZero();
  u(3) = 9.81*m;
  vector<VectorXd> us(N, u);

  AirmDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = 1;

  AirmView view(sys, &ddp.xs);
  if (viewer)
    viewer->Add(view);  

  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed

  getchar();

  for (int i = 0; i < 100; ++i) {    
    timer_start(timer);
    ddp.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << ": took " << te << " us." << endl;        
    getchar();
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
  viewer->frameName = "../../logs/airm/frames/frame";

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
