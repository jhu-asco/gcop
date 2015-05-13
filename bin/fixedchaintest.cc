#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "chainview.h"
#include "utils.h"
#include "lqcost.h"
#include "params.h"
#include "mbscontroller.h"
#include "se3.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<MbsState> ChainDdp;

Params params;

void solver_process(Viewer* viewer)
{

  if (viewer) {
    viewer->SetCamera(-25.8, 56, -0.14, 0, -0.63);
  }

  int N = 128;      // discrete trajectory segments
  double tf = 5;   // time-horizon

  params.GetInt("N", N);
  params.GetDouble("tf", tf);

  int nb = 4;     // nof bodies
  params.GetInt("nb", nb);

  double h = tf/N; // time-step

  // fixed chain
  Chain sys(nb,true);
  sys.debug=false;

  params.GetInt("method", sys.method);
  params.GetInt("iters", sys.iters);

  params.GetVectorXd("damping", sys.damping);
  
  params.GetVectorXd("ulb", sys.U.lb);
  params.GetVectorXd("uub", sys.U.ub);

  // acceleration due to gravity
  params.GetVector3d("ag", sys.ag);
	sys.Init();

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

  LqCost<MbsState> cost(sys, tf, xf);
  
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

  // not using it @MK: this is the new part, initialize trajectory using a controller
  MbsController ctrl(sys);
  params.GetVectorXd("Kp", ctrl.Kp);
  params.GetVectorXd("Kd", ctrl.Kd);

  for (int i = 0; i < xs.size()-1; ++i) {
    double t = i*h;
    // ctrl.Set(us[i], t, xs[i]); 
    sys.Step(xs[i+1], i*h, xs[i], us[i], h);
  }

  // see the result before running optimization
  //  getchar();

  ChainDdp ddp(sys, cost, ts, xs, us);
  params.GetDouble("mu", ddp.mu);

  params.GetDouble("eps", ddp.eps);

  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed  

  for (int i = 0; i < 50; ++i) {
    timer_start(timer);
    ddp.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << ": took " << te << " us." << endl;
    //    getchar();
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
    params.Load("../../bin/fixedchain1.cfg");


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
  float rgb[3] = {1,1,1};
  viewer->SetColor(rgb);
  viewer->Start();
#endif


  return 0;
}

