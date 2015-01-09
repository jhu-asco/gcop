#include <iomanip>
#include <iostream>
#include "gndocp.h"
#include "viewer.h"
#include "rccarview.h"
#include "utils.h"
#include "rnlqcost.h"
#include "params.h"

#include "controltparam.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef GnDocp<Vector4d, 4, 2, Dynamic, 6> RccarGn;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(-5, 51, -0.2, -0.15, -2.3);

  int N = 32;        // number of segments
  double tf = 5;    // time horizon

  int iters = 30;

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);

  params.GetInt("iters", iters);
  

  double h = tf/N;   // time step

  Rccar sys;

  //  sys.U.lb[1] = tan(-M_PI/5);
  //  sys.U.ub[1] = tan(M_PI/5);

  // initial state
  Vector4d x0(1,1,0,0);
  params.GetVector4d("x0", x0);

  // final state
  Vector4d xf(0,0,0,0);
  params.GetVector4d("xf", xf);  

  // cost
  RnLqCost<4, 2, Dynamic, 6> cost(sys, tf, xf);
  VectorXd Q(4);
  if (params.GetVectorXd("Q", Q))
    cost.Q = Q.asDiagonal();

  VectorXd Qf(4);
  if (params.GetVectorXd("Qf", Qf))
    cost.Qf = Qf.asDiagonal();
  
  VectorXd R(2);
  if (params.GetVectorXd("R", R)) 
    cost.R = R.asDiagonal();

  cost.UpdateGains();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Vector4d> xs(N+1);
  // initial state
  xs[0] = x0;

  // initial controls
  vector<Vector2d> us(N);
  for (int i = 0; i < N/2; ++i) {
    us[i] = Vector2d(.01, .0);
    us[N/2+i] = Vector2d(-.01, .0);
  }
  
  //  Tparam<Vector4d, 4, 2> tp(sys, us.size()*sys.U.n);

  int Nk = 10;
  vector<double> tks(Nk+1);
  for (int k = 0; k <=Nk; ++k)
    tks[k] = k*(tf/Nk);
  
  ControlTparam<Vector4d, 4, 2> ctp(sys, tks);

  RccarGn gn(sys, cost, ctp, ts, xs, us);  

  RccarView view(sys, &gn.xs);
  
  viewer->Add(view);

  struct timeval timer;
  // Gn.debug = false; // turn off debug for speed
  getchar();

  for (int i = 0; i < iters; ++i) {
    timer_start(timer);
    gn.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << " took: " << te << " us." << endl;
    cout<<"Optimal Cost: "<<gn.J<<endl;
    getchar();
  }

  cout << xs[N] << endl;

  //  xs[1][3]  velocity
  //atan(us[0][1]) steering angle
 
  cout << "done!" << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{

  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("../../bin/rccar.cfg"); 


#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "../../logs/rccar/frames/frame";

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
