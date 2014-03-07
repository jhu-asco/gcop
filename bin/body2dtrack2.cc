#include <iomanip>
#include <iostream>
#include "body2dslam.h"
#include "utils.h"
#include "se2.h"
#include "viewer.h"
#include "body2dview.h"
#include "body2dcost.h"
#include "params.h"
#include "body2dtrackview.h"

using namespace std;
using namespace Eigen;
using namespace gcop;


typedef Dmoc<pair<Matrix3d, Vector3d>, 6, 3> Body2dDmoc;

Params params;

void Run(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(18.875, 1.625, -4, -1.75, -5.45);

  int N = 100;
  double tf = 10;
  int nf = 200;
  params.GetInt("N", N);  
  params.GetDouble("tf", tf);
  params.GetInt("nf", nf);
  double h = tf/N;

  // control horizon params
  int iters = 30;
  double Tc = 2;
  params.GetInt("iters", iters);
  params.GetDouble("Tc", Tc);
  int Nc = (int)(Tc/h);

  Body2d sys(new Body2dForce(true));
  sys.force->D << .01, .01, 3; // add damping

  // options: odometry, paramForce, forces
  Body2dTrack pgt(sys, N, nf, 0, tf, false, false, true);   ///< ground truth
  pgt.MakeTrue();

  Body2dTrack pg(sys, N, nf, 0, tf, false, false, true);    ///< noisy pose track

  //  Body2dTrack::Synthesize3(pgt, pg, tf);

  Body2dView tview(sys, &pgt.xs);   // true path
  tview.lineWidth = 5;
  tview.rgba[0] = 0;
  tview.rgba[1] = 1;
  tview.rgba[2] = 0;
  tview.renderSystem = false;

  Body2dTrackView pgtv(pgt);             // true pose-track
  pgtv.rgba[0] = 0;   pgtv.rgba[1] = 1;
  //  pgtv.drawLandmarks = false;

  Body2dView view(sys, &pg.xs);     // optimized path
  view.lineWidth = 5;
  view.rgba[0] = 1;
  view.rgba[1] = 0;
  view.rgba[2] = 0;
  view.renderSystem = false;


  Body2dTrackView pgv(pg);               // optimized pose-track
  pgv.rgba[0] = 1; pgv.rgba[0] = 0;
  pgv.drawLandmarks = false;

  if (viewer) {
    viewer->Add(tview);
    viewer->Add(pgtv);
    
    viewer->Add(view);
    viewer->Add(pgv);
  }

  struct timeval timer;
  getchar();
  

  SE2 &se2 = SE2::Instance();  

  M3V3d x0 = pgt.xs[0];  
  M3V3d xf = pgt.xs[Nc];
  cout << xf.first << endl;
  cout << xf.second << endl;

  // cost
  Body2dCost cost(Tc, xf);

  VectorXd Q(6);
  if (params.GetVectorXd("Q", Q))
    cost.Q = Q.asDiagonal();
  
  VectorXd Qf(6);
  if (params.GetVectorXd("Qf", Qf))
    cost.Qf = Qf.asDiagonal();
  
  VectorXd R(3);
  if (params.GetVectorXd("R", R)) 
    cost.R = R.asDiagonal();

  // times
  vector<double> ts(Nc+1);
  for (int k = 0; k <=Nc; ++k)
    ts[k] = k*h;

  // states
  vector<pair<Matrix3d, Vector3d> > xs(Nc+1, x0);

  // controls
  Vector3d u(0,0,0);
  vector<Vector3d> us(Nc, u);

  // past trajectory
  vector<double> tps(Nc+1);
  for (int k = 0; k <=Nc; ++k)
    tps[k] = -k*h;
  vector<pair<Matrix3d, Vector3d> > xps(Nc+1, x0);
  vector<Vector3d> ups(Nc, u);


  Body2dDmoc dmoc(sys, cost, ts, xs, us);
  dmoc.mu = .01;
  params.GetDouble("mu", dmoc.mu);

  Body2dView cview(sys, &dmoc.xs);
  viewer->Add(cview);

  Body2dView pview(sys, &xps);
  viewer->Add(pview);
  pview.rgba[0] = 1;  pview.rgba[1] = 1;  pview.rgba[2] = 0;
  
  for (int i=Nc; i<=N; ++i) {
    
    xf = pgt.xs[i];

    for (int j = 0; j < iters; ++j) {      
      timer_start(timer);
      dmoc.Iterate();
      long te = timer_us(timer);      
      cout << "Iteration #" << j << " took: " << te << " us." << endl;    
    }

    // shift control trajectory forward
    for (int k = 0; k < Nc; ++k) {
      ts[k] = ts[k+1];
      if (k < Nc - 1)
        us[k] = us[k+1];

      tps[k] = tps[k+1];
      xps[k] = xps[k+1];
      if (k < Nc - 1)
        ups[k] = ups[k+1];

    }
    ts[Nc] = ts[Nc] + h;
    xs[0] = xs[1];

    tps[Nc] = ts[0];
    xps[Nc] = xs[0];


    getchar();  
  }


  /*
  Body2dSlam ba(pg);
  ba.pdmoc->debug = true; // turn off debug for speed
  ba.pdmoc->mu = .01;
  ba.pdmoc->nu = .01;

  for (int i = 0; i < 1000; ++i) {

    cout << "Press Enter to continue" << endl;
    getchar();    
    
    timer_start(timer);
    ba.pdmoc->Iterate();
    long te = timer_us(timer);
    cout << ba.pdmoc->dus[0] << endl;
    cout << "Iteration #" << i << " took: " << te << " us." << endl;    
    cout << "p=" << ba.pdmoc->p.head<2>().transpose() << endl;    
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
    params.Load("../../bin/body2dtrack.cfg");

#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "videos/sys";

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
