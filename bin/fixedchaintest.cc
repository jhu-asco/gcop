#include <iomanip>
#include <iostream>
#include "dmoc.h"
#include "viewer.h"
#include "chainview.h"
#include "utils.h"
#include "lqcost.h"
#include "params.h"
#include "mbscontroller.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<MbsState> ChainDmoc;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer) {
    viewer->SetCamera(45.25, 37, -0.14, 0.05, -1.75);
  }

  int N = 128;      // discrete trajectory segments
  double tf = 5;   // time-horizon

  params.GetInt("N", N);
  params.GetDouble("tf", tf);  

  int nb = 3;     // nof bodies
  params.GetInt("nb", nb);

  double h = tf/N; // time-step

  // fixed chain
  Chain sys(nb, true);

  params.GetInt("method", sys.method);
  params.GetInt("iters", sys.iters);

  params.GetVectorXd("damping", sys.damping);
  
  params.GetVectorXd("ulb", sys.U.lb);
  params.GetVectorXd("uub", sys.U.ub);

  // acceleration due to gravity
  params.GetVector3d("ag", sys.ag);

  int n = nb - 1; // movable bodies

  VectorXd qv0(2*n);
  params.GetVectorXd("x0", qv0);

  MbsState x0(nb, true);
  x0.gs[0].setIdentity();    
  x0.r = qv0.head(n);
  x0.dr = qv0.tail(n);
  sys.Rec(x0, h);

  VectorXd qvf(2*n);
  params.GetVectorXd("xf", qvf);  

  MbsState xf(nb, true);
  xf.gs[0].setIdentity();
    
  xf.r = qvf.head(n);
  xf.dr = qvf.tail(n);

  LqCost<MbsState> cost(sys.X, (Rn<>&)sys.U, tf, xf);
  
  VectorXd Q(2*n);
  VectorXd R(n);
  VectorXd Qf(2*n);

  params.GetVectorXd("Q", Q);
  params.GetVectorXd("R", R);
  params.GetVectorXd("Qf", Qf);
 
  cost.Q = Q.asDiagonal();
  cost.R = R.asDiagonal();
  cost.Qf = Qf.asDiagonal();

  // times
  vector<double> ts(N + 1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // states
  vector<MbsState> xs(N + 1, x0);

  // initial controls
  VectorXd u(n);
  u.setZero();
  vector<VectorXd> us(N, u);

  ChainView view(sys, &xs);
  if (viewer)
    viewer->Add(view);

  // @MK: this is the new part, initialize trajectory using a controller
  MbsController ctrl(sys);
  params.GetVectorXd("Kp", ctrl.Kp);
  params.GetVectorXd("Kd", ctrl.Kd);

  for (int i = 0; i < xs.size()-1; ++i) {
    double t = i*h;
    // ctrl.Set(us[i], t, xs[i]); 
    sys.Step(xs[i+1], i*h, xs[i], us[i], h);
  }

  // see the result before running optimization
  getchar();

  ChainDmoc dmoc(sys, cost, ts, xs, us);
  params.GetDouble("mu", dmoc.mu);

  params.GetDouble("eps", dmoc.eps);

  struct timeval timer;
  //  dmoc.debug = false; // turn off debug for speed

  for (int i = 0; i < 50; ++i) {
    timer_start(timer);
    dmoc.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << ": took " << te << " us." << endl;
    getchar();
  }

  int Nd;
  params.GetInt("Nd", Nd);
  vector<MbsState> xds(Nd+1, x0);
  int d = N/Nd;
  for (int i=Nd, j = N; i>=0 && j>= 0; --i, j-=d)
    xds[i] = xs[j];

  //  xs = xds;
  
  cout << "done!" << endl;
  cout << "Click on window and press 'a' to animate." << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{

  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("../../bin/fixedchain.cfg");


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
  viewer->SetColor((float[3]){1,1,1});
  viewer->Start();
#endif


  return 0;
}

