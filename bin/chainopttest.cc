#include <iomanip>
#include <iostream>
#include "dmoc.h"
#include "viewer.h"
#include "chainview.h"
#include "utils.h"
#include "lqcost.h"
#include "rn.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<MbsState> ChainDmoc;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(45.25, 37, -0.14, 0.05, -1.75);

  int N = 64;      // discrete trajectory segments
  double tf = 20;   // time-horizon
  double h = tf/N; // time-step

  // system
  Chain sys;

  MbsState xf(3);
  xf.gs[0].setIdentity();
  xf.r[0] = M_PI/4;   
  xf.r[1] = M_PI/4;
  xf.vs[0].setZero();
  xf.dr.setZero();

  LqCost<MbsState> cost(sys.X, (Rn<>&)sys.U, tf, xf);

  cost.Qf(0,0) = 2; cost.Qf(1,1) = 2; cost.Qf(2,2) = 2;
  cost.Qf(3,3) = 20; cost.Qf(4,4) = 20; cost.Qf(5,5) = 20;
  //  cost.R *= .0;
  //  cost.Q.setZero();


  /*
  cost.Qf(0,0) = 2; cost.Qf(1,1) = 2; cost.Qf(2,2) = 2;
  cost.Qf(3,3) = 20; cost.Qf(4,4) = 20; cost.Qf(5,5) = 20;
  cost.Qf(6,6) = 1; cost.Qf(7,7) = 1; 

  cost.Qf(8,8) = 5; cost.Qf(9,9) = 5; cost.Qf(10,10) = 5;
  cost.Qf(11,11) = 50; cost.Qf(12,12) = 50; cost.Qf(13,13) = 50;
  cost.Qf(14,14) = 5; cost.Qf(15,15) = 5; 

  cost.R(0,0) = .05; cost.R(1,1) = .05; cost.R(2,2) = .05; 
  cost.R(3,3) = .1; cost.R(4,4) = .1; cost.R(5,5) = .1;
  cost.R(6,6) = .01;  cost.R(7,7) = .01;

  */

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  MbsState x(3);
  SE3::Instance().rpyxyz2g(x.gs[0], Vector3d(M_PI/4,M_PI/4,0), Vector3d(-2,0,0));
  x.r[0] = -M_PI/2;   
  x.r[1] = -M_PI/2;
  x.vs[0].setZero();
  x.dr[0] = 0;
  x.dr[1] = 0;
  sys.Rec(x, h);

  // states
  vector<MbsState> xs(N+1, x);
  //  xs[0].vs[0].setZero();
  //  xs[0].dr[0] = 0;
  //  xs[0].dr[1] = 0;
  //  sys.KStep(xs[0], x, h);


  double m = sys.links[0].m + sys.links[1].m  + sys.links[2].m;
  
  // initial controls (e.g. hover at one place)
  VectorXd u(8);
  u.setZero();
	u[5] = m*9.81;
  vector<VectorXd> us(N, u);

  ChainDmoc dmoc(sys, cost, ts, xs, us);
  dmoc.mu = 10;

  ChainView view(sys, &dmoc.xs);
  if (viewer)
    viewer->Add(view);  



  struct timeval timer;
  //  dmoc.debug = false; // turn off debug for speed

  for (int i = 0; i < 100; ++i) {    
    timer_start(timer);
    dmoc.Iterate();
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
  viewer->frameName = "../../logs/chain/frames/frame";

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
