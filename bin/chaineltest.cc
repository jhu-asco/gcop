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

  double h = tf/N; // time-step

  // system
  Chain sys;

  params.GetInt("method", sys.method);  
  params.GetInt("iters", sys.iters);  
 

  // acceleration due to gravity
  params.GetVector3d("ag", sys.ag);

  VectorXd qv0(16);
  params.GetVectorXd("x0", qv0);  

  MbsState x0(3);
  SE3::Instance().rpyxyz2g(x0.gs[0], qv0.head(3), qv0.segment<3>(3));
    
  x0.r = qv0.segment<2>(6);
  x0.vs[0] = qv0.segment<6>(8);
  x0.dr = qv0.tail(2);
  sys.Rec(x0, h);

  VectorXd qvf(16);
  params.GetVectorXd("xf", qvf);  

  MbsState xf(7);
  SE3::Instance().rpyxyz2g(xf.gs[0], qvf.head(3), qvf.segment<3>(3));
    
  xf.r = qvf.segment<2>(6);
  xf.vs[0] = qvf.segment<6>(8);
  xf.dr = qvf.tail(2);

  LqCost<MbsState> cost(sys.X, (Rn<>&)sys.U, tf, xf);
  
  VectorXd Q(16);
  VectorXd R(8);
  VectorXd Qf(16);

  params.GetVectorXd("Q", Q);
  params.GetVectorXd("R", R);
  params.GetVectorXd("Qf", Qf);
 
  cost.Q = Q.asDiagonal();
  cost.R = R.asDiagonal();
  cost.Qf = Qf.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // states
  vector<MbsState> xs(N+1, x0);

  double m = sys.links[0].m + sys.links[1].m  + sys.links[2].m;
	cout<<"Mass: "<<m<<endl;
  
  // initial controls (e.g. hover at one place)
  VectorXd u(8);
  u.setZero();
  //  u(5) = 9.81*m;
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
    //    ctrl.Set(us[i], t, xs[i]); 
    sys.Step(xs[i+1], i*h, xs[i], us[i], h);
  }

  // see the result before running optimization
  getchar();

  ChainDmoc dmoc(sys, cost, ts, xs, us);
  params.GetDouble("mu", dmoc.mu);

  struct timeval timer;
  dmoc.debug = false; // turn off debug for speed

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
    params.Load("../../bin/chaineltest.cfg");


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

