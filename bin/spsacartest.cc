#include <iomanip>
#include <iostream>
#include "spsa.h"
#include "viewer.h"
#include "rccarview.h"
#include "utils.h"
#include "rnlqcost.h"
#include "params.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef SPSA<Vector4d, 4, 2> RccarSPSA;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(-5, 51, -0.2, -0.15, -2.3);

  int N = 10;        // number of segments
  double tf = 5;    // time horizon

  int iters = 30;

  params.GetInt("N", N);  

  params.GetInt("iters", iters);  

  params.GetDouble("tf", tf);


  double h = tf/N;   // time step

  Rccar sys;

  //  sys.U.lb[1] = tan(-M_PI/5);
  //  sys.U.ub[1] = tan(M_PI/5);

  Vector4d stepcfs(0,0,0,0);
  params.GetVector4d("stepcoeffs", stepcfs); //Step coeffs for SPSA

  // initial state
  Vector4d x0(1,1,0,0);
  params.GetVector4d("x0", x0);

  // final state
  Vector4d xf(0,0,0,0);
  params.GetVector4d("xf", xf);  

  // cost
  RnLqCost<4, 2> cost(tf, xf);
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
	Vector2d u0;
	params.GetVector2d("u0",u0);

  for (int i = 0; i < N/2; ++i) {
    us[i] = Vector2d(u0(0), u0(1));
    us[N/2+i] = Vector2d(-u0(0), -u0(1));    
  }

  RccarSPSA spsa(sys, cost, ts, xs, us);
	spsa.debug = false;
  //  dmoc.mu = .01;
  //  params.GetDouble("mu", dmoc.mu);

  params.GetInt("Nit", spsa.Nit);

	spsa.stepc.a = stepcfs[0];
	spsa.stepc.c1 = stepcfs[1];
	spsa.stepc.alpha = stepcfs[2];
	spsa.stepc.gamma = stepcfs[3];//Set coeffs based on parameters
	spsa.stepc.A = 0.1*spsa.Nit*iters;//10 percent of total number of iterations


  RccarView view(sys, &spsa.xs);
  
  viewer->Add(view);

  struct timeval timer;
  // dmoc.debug = false; // turn off debug for speed
  getchar();

  for (int i = 0; i < iters; ++i) {
    timer_start(timer);
    spsa.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << " took: " << te << " us." << endl;
    cout << "Cost=" << spsa.J << endl;
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
    params.Load("../../bin/spsacar.cfg");


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
