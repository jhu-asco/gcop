#include <iomanip>
#include <iostream>
#include "systemce.h"
#include "viewer.h"
#include "rccarview.h"
#include "utils.h"
#include "rnlqcost.h"
#include "params.h"

#define USE_SPLINEPARAM

#ifdef USE_SPLINEPARAM
#include "uniformsplinetparam.h"
#else
#include "controltparam.h"
#endif

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef SystemCe<Vector4d, 4, 2, Dynamic> RccarCe;

Params params;

// whether to use a specific trajectory parametrization
#define USE_TPARAM

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

  cout<<"sys.U.n: "<<sys.U.n<<endl;

  //  sys.U.lb[1] = tan(-M_PI/5);
  //  sys.U.ub[1] = tan(M_PI/5);

  // initial state
  Vector4d x0(1,1,0,0);
  params.GetVector4d("x0", x0);

  // final state
  Vector4d xf(0,0,0,0);
  params.GetVector4d("xf", xf);  

  // cost
  RnLqCost<4, 2> cost(sys, tf, xf);
  VectorXd Q(4);
  if (params.GetVectorXd("Q", Q))
    cost.Q = Q.asDiagonal();
  
  VectorXd Qf(4);
  if (params.GetVectorXd("Qf", Qf))
    cost.Qf = Qf.asDiagonal();
  
  VectorXd R(2);
  if (params.GetVectorXd("R", R)) 
    cost.R = R.asDiagonal();

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

  Vector2d du(.5, .33);
  params.GetVector2d("du", du);

  Vector2d e(.001, .001);
  params.GetVector2d("e", e);

  vector<Vector2d> dus(N, du);
  vector<Vector2d> es(N, e);
  
#ifdef USE_TPARAM
  int Nk = 5;
  //vector<double> tks(Nk+1);
  VectorXd tks(Nk+1);
  for (int k = 0; k <=Nk; ++k)
  {
    double u = double(k)/Nk;
    //tks[k] = tf*u*u;
    tks[k] = tf*u;
  }
  

#ifdef USE_SPLINEPARAM
  UniformSplineTparam<Vector4d, 4, 2> ctp(sys, tks);
#else
  ControlTparam<Vector4d, 4, 2> ctp(sys, tks);
#endif

  RccarCe ce(sys, cost, ctp, ts, xs, us, 0, dus, es);
#else
  RccarCe ce(sys, cost, ts, xs, us, 0, dus, es);
#endif

  params.GetBool("mras", ce.ce.mras);
  params.GetBool("inc", ce.ce.inc);

  //  ddp.mu = .01;
  //  params.GetDouble("mu", ddp.mu);

  params.GetInt("Ns", ce.Ns);


  RccarView view(sys, &ce.xs);
  
  viewer->Add(view);

  struct timeval timer;
  // ddp.debug = false; // turn off debug for speed
  getchar();

  for (int i = 0; i < iters; ++i) {
    timer_start(timer);
    ce.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << " took: " << te << " us." << endl;
    cout << "Cost=" << ce.J << endl;
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
    params.Load("../../bin/cecar.cfg");


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
