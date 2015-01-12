#include <iomanip>
#include <iostream>
#include "body3dtrackcost.h"
#include "pddp.h"
#include "utils.h"
#include "se2.h"
#include "viewer.h"
#include "body3dview.h"
#include "body3dcost.h"
#include "params.h"
#include "body3dtrackview.h"

using namespace std;
using namespace Eigen;
using namespace gcop;


typedef Ddp<Body3dState, 12, 6> Body3dDdp;

Params params;

void Run(Viewer* viewer)
{
  srand(1);
  
  if (viewer)
    viewer->SetCamera(18.875, 1.625, -.15, -.6, -35.5);
  
  double tf = 30;
  int nf = 100;
  params.GetDouble("tf", tf);
  params.GetInt("nf", nf);

  double Ts = 1;
  params.GetDouble("Ts", Ts);

  bool renderForces = false;
  params.GetBool("renderForces", renderForces);

  bool hideTrue = false;
  params.GetBool("hideTrue", hideTrue);
  bool hideEst = false;
  params.GetBool("hideEst", hideEst);

  bool hideOdom = false;
  params.GetBool("hideOdom", hideOdom);


  // control horizon params
  int iters = 30;
  double Tc = 2; // horizon
  int N; // optimal control segments
  params.GetInt("iters", iters);
  params.GetDouble("Tc", Tc);
  params.GetInt("N", N);
  double h = Tc/N;

  Body3d<> sys;

  double r = 25;
  params.GetDouble("r", r);

  double vd = 5;
  params.GetDouble("vd", vd);

  // options: odometry, paramForce, forces
  Body3dTrack pg(sys, nf, 0, tf, r, true, false, true);   ///< ground truth

  params.GetDouble("w", pg.w);
  params.GetVector6d("cw", pg.cw);
  params.GetVector6d("cv", pg.cv);
  params.GetDouble("cp", pg.cp);
  params.GetDouble("dmax", pg.dmax);

  pg.MakeTrue();

  Body3dState x0;
  pg.Get(x0, vd, 0);

  //  Body2dTrack pg(sys, N, nf, 0, tf, false, false, true);    ///< noisy pose track
  //  Body2dTrack::Synthesize3(pgt, pg, tf);

  Body3dView<> tview(sys, &pg.xs, &pg.us);   // estimated path
  tview.lineWidth = 5;
  tview.rgba[0] = 0;
  tview.rgba[1] = 1;
  tview.rgba[2] = 0;
  tview.renderSystem = renderForces;
  tview.renderForces = renderForces;

  Body3dView<> oview(sys, &pg.xos);   // odometry
  oview.lineWidth = 5;
  oview.renderSystem = false;
  oview.rgba[0] = 1;
  oview.rgba[1] = 0;
  oview.rgba[2] = 0;

  Body3dTrackView pgv(pg);               // optimized pose-track
  pgv.rgba[0] = 1; pgv.rgba[1] = 1; pgv.rgba[2] = 1; 

  pgv.drawLandmarks = true;

  pgv.drawForces = false;
  pgv.forceScale=.1;


  if (viewer) {
    if (!hideEst)
      viewer->Add(tview);
    if (!hideOdom)
      viewer->Add(oview);
    viewer->Add(pgv);
  }

  struct timeval timer;
  getchar();
  

  SE2 &se2 = SE2::Instance();  
  Body3dState xf;
  pg.Get(xf, vd, Tc);

  // cost
  Body3dCost<> cost(sys, Tc, xf);
  //params.GetDouble("ko", cost.ko);
  //if (cost.ko > 1e-10)
  //  cost.track = &pg;


  VectorXd Q(12);
  if (params.GetVectorXd("Q", Q))
    cost.Q = Q.asDiagonal();
  
  VectorXd Qf(12);
  if (params.GetVectorXd("Qf", Qf))
    cost.Qf = Qf.asDiagonal();
  
  VectorXd R(6);
  if (params.GetVectorXd("R", R)) 
    cost.R = R.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Body3dState > xs(N+1, pg.xs[0]);

  // controls
  Vector6d u;
  u << 0,0,0,0,0,0;
  vector<Vector6d> us(N, u);

  // past true trajectory  
  vector<double> tps(1,0);
  vector<Body3dState > xps(1, x0);
  vector<Vector6d> ups;

  Body3dDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = .01;
  params.GetDouble("mu", ddp.mu);

  Body3dView<> cview(sys, &ddp.xs);
  cview.rgba[0] = 0;  cview.rgba[1] = 1;  cview.rgba[2] = 1;
  viewer->Add(cview);
  cview.renderSystem = false;
  cview.lineWidth = 4;

  // cview.renderForces = true;

  Body3dView<> pview(sys, &xps, &ups);
  if (!hideTrue)
      viewer->Add(pview);
  pview.rgba[0] = 1;  pview.rgba[1] = 1;  pview.rgba[2] = 0;
  pview.renderSystem = renderForces;
  pview.renderForces = renderForces;

  Body3dTrackCost tcost(0, pg);  ///< cost function
  PDdp<Body3dState, 12, 6> *pddp = 0;


  for (double t=0; t < tf; t+=h) {

    pg.Get(xf, vd, t+Tc);
    cost.tf = t + Tc;

    bool oc = false;
    params.GetBool("oc", oc);

    xs[0] = pg.xs.back();    
    for (int j=0; j < N; ++j) {
      sys.Step(xs[j+1], t + j*h, xs[j], us[j], h);
    }

    if (oc) {
      for (int j = 0; j < iters; ++j) {
        timer_start(timer);
        ddp.Iterate();
        long te = timer_us(timer);      
        cout << "Iteration #" << j << " took: " << te << " us." << endl;    
      }
    }

    // apply actual control (puluated with noise)
    Vector6d w;
    w << sqrt(pg.cw[0])*random_normal(), 
                          sqrt(pg.cw[1])*random_normal(), 
                          sqrt(pg.cw[2])*random_normal(), 
                          sqrt(pg.cw[3])*random_normal(), 
                          sqrt(pg.cw[4])*random_normal(), 
                          sqrt(pg.cw[5])*random_normal(); 

    //w << 0, sqrt(pg.cw[1]), 0, 0, 0, 0;
    
    // simulate true state
    Body3dState xt;
    sys.Step(xt, t, xps.back(), us[0] + w, h);

    // add assumed control and true state to estimator
    pg.Add2(us[0], xt, h);

    // shift control trajectory forward
    for (int k = 0; k < N; ++k) {
      ts[k] = ts[k+1];
      if (k < N - 1)
        us[k] = us[k+1];
    }
    ts[N] = ts[N] + h;

    tcost.tf = t+h;
    
    if (t > Ts) {
      cout << "ts " << pg.ts.size() << endl;
      cout << "xs " << pg.xs.size() << endl;
      cout << "us " << pg.us.size() << endl;
      cout << "p " << pg.p.size() << endl;

      pddp = new PDdp<Body3dState, 12, 6>(pg.sys, tcost, pg.ts, pg.xs, pg.us, pg.p, 3*pg.extforce);
      for (int b=0; b < 20;++b)
        pddp->Iterate();     

      delete pddp;
    }

    tps.push_back(ts[0]);
    xps.push_back(xt);
    ups.push_back(us[0] + w);
    //    ups.push_back(us[0] - pg.us[0]);

    cout << "FEATURES:" << pg.p.size() << endl;
    
    //getchar(); 
    //viewer->saveSnapshot=true;
  }


  /*
  Body2dSlam ba(pg);
  ba.pddp->debug = true; // turn off debug for speed
  ba.pddp->mu = .01;
  ba.pddp->nu = .01;

  for (int i = 0; i < 1000; ++i) {

    cout << "Press Enter to continue" << endl;
    getchar();    
    
    timer_start(timer);
    ba.pddp->Iterate();
    long te = timer_us(timer);
    cout << ba.pddp->dus[0] << endl;
    cout << "Iteration #" << i << " took: " << te << " us." << endl;    
    cout << "p=" << ba.pddp->p.head<2>().transpose() << endl;    
    cout << "ut=" << pgt.us[2] << endl;    
    cout << "u=" << pg.us[2] << endl;    

  }
  */
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
    params.Load("../../bin/body3dtrack.cfg");

#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "../../logs/body3dtrack/frames/body3d";
  viewer->displayName = "../../logs/body3dtrack/display/body3d";

  pthread_t dummy;
  pthread_create( &dummy, NULL, (void *(*) (void *)) Run, viewer);

#else
  Run(0);
#endif


#ifdef DISP
  viewer->Start();
#endif


  return 0;
}
