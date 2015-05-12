#include <iostream>
#include "viewer.h"
#include "body3dview.h"
#include "body3davoidcontroller.h"
#include "utils.h"
#include "disk.h"
#include "systemce.h"
#include "params.h"
#include "lqcost.h"
#include "controllertparam.h"
#include "demview.h"
#include "pqpdem.h"
#include "systemceview.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

//typedef Ddp<Body3dState, 4, 2> Body3dDdp;
//typedef Disk<Body3dState, 4, 2> Body3dDisk;
//typedef ConstraintCost<Body3dState, 4, 2, Dynamic, 1> DiskCost;

typedef SystemCe<Body3dState, 12, 6, Dynamic, 5> Body3dCe;
typedef SystemCeView<Body3dState, 12, 6, Dynamic, 5> Body3dCeView;

typedef PqpDem<Body3dState, 12, 6> Body3dPqpDem;

Params params;

class Body3dSampler : public Sampler<Body3dState> {
public:
  Body3dSampler(Body3dState x0, double r, int N): i(0) {
    Body3dState x;
    for (int i = 0; i < N; ++i) {
      
      x.second[0] = 180;
      x.second[1] = 180*RND;
      x.second[2] = 2;      

      /*
      double a = RND*M_PI - M_PI/2;
      x.second[0] = x0.second[0] + r*cos(a);
      x.second[1] = x0.second[1] + r*sin(a);
      x.second[2] = x0.second[2];      
      */
      samples.push_back(x);
    }
  }

  bool operator()(Body3dState &x) {
    x = samples[i%samples.size()];
    i++;
    return true;
  }

  int i;
  vector<Body3dState> samples;
};

void solver_process(Viewer* viewer)
{
  viewer->SetCamera(-7.25, 61, -2.8, -2.55, -10);

  Body3d<6> sys;

  int N = 256;
  double tf = 5;
  int iters = 100;

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);
  params.GetInt("iters", iters);

  // cost 
  Body3dState x0(Matrix3d::Identity(), Vector9d::Zero());
  x0.second[0] = 46;
  x0.second[1] = 82;
  x0.second[2] = 2;

  Body3dState xf(Matrix3d::Identity(), Vector9d::Zero());
  xf.second[0] = 160;
  xf.second[1] = 125;  
  xf.second[2] = 2;

  VectorXd qv0(12);
  params.GetVectorXd("x0", qv0);  
  SO3::Instance().q2g(x0.first, qv0.head(3));    
  x0.second = qv0.segment<9>(3);

  vector<Body3dState> xs(N+1);
  xs[0] = x0;

  VectorXd qvf(12);
  params.GetVectorXd("xf", qvf);  
  SO3::Instance().q2g(xf.first, qvf.head(3));    
  xf.second = qvf.segment<9>(3);
  
  float camParams[5];
  if (viewer) {
    params.GetFloatArray("camParams", 5, camParams);
    viewer->SetCamera(camParams);
  }

  string mapName("maps/pic2_dilated.ppm");
  double mapCellSize = .25;
  double mapHeightScale = 30;
  params.GetString("mapName", mapName);
  params.GetDouble("mapCellSize", mapCellSize);
  params.GetDouble("mapHeightScale", mapHeightScale);
  Dem dem(mapName.c_str(), mapCellSize, mapHeightScale);
  //  dem.Convolve(1);

  double cr = 0;
  params.GetDouble("cr", cr);

  Body3dPqpDem pqp(dem, cr);  

  Body3dAvoidController<> ctrl(sys, &xf, 0, &pqp);
  ctrl.avoidCtrl->k = 5; // avoidance gain
  ctrl.avoidCtrl->sr = 20; // avoidance gain

  params.GetDouble("kb", ctrl.avoidCtrl->kb);
  params.GetDouble("sr", ctrl.avoidCtrl->sr);

  DemView dv(dem);
  if (viewer)
    viewer->Add(dv);

  // goal view
  vector<Body3dState> xfs;
  xfs.push_back(xf);
  Body3dView<6> goalView(sys, &xfs);
  goalView.rgba[0] = 1;
  if (viewer) 
    viewer->Add(goalView);

  LqCost<Body3dState, 12, 6> cost(sys, tf, xf);
  VectorXd Q(12);
  VectorXd R(6);  
  VectorXd Qf(12);
  params.GetVectorXd("Q", Q);
  params.GetVectorXd("R", R);
  params.GetVectorXd("Qf", Qf);
  cost.Q = Q.asDiagonal();
  cost.R = R.asDiagonal();
  cost.Qf = Qf.asDiagonal();


  // times
  double h = tf/N;  
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;
  
  // initial controls
  vector<Vector6d> us(N);
  
  ControllerTparam<Body3dState, 12, 6, Dynamic, 5, Body3dState> ctp(sys, ctrl, 5, &pqp);
  
  Vector5d mu0;
  params.GetVector5d("mu0", mu0); 

  Vector5d P0d;
  params.GetVector5d("P0", P0d);
  Matrix5d P0 = P0d.asDiagonal();

  Vector5d Sd;
  params.GetVector5d("S", Sd);
  Matrix5d S = Sd.asDiagonal();

  // process noise terms
  params.GetDouble("sw", sys.sw);
  params.GetDouble("sv", sys.sv);

  Body3dCe ce(sys, cost, ctp, ts, xs, us, 0, mu0, P0, S);
  ce.ce.gmm.ns[0].bounded = true;
  ce.ce.gmm.ns[0].lb.setZero();
  ce.ce.gmm.ns[0].ub.setConstant(1000);
  ce.Ns = 200;
  params.GetInt("Ns", ce.Ns);

  ctp.contextSampler = new Body3dSampler(x0, 100, ce.Ns);
  
  ce.Jub = 1000000;
  ce.enforceUpperBound = true;

  Body3dView<> view(sys, &xs);
  viewer->Add(view);  
  view.rgba[0] = 0; view.rgba[1] = 0; view.rgba[2] = 1;
  
  Body3dCeView ceview(ce, view);
  viewer->Add(ceview);

  struct timeval timer;

  getchar();

  for (int i = 0; i < 100; ++i) {
    timer_start(timer);
    //    rseed(0);
    srand(0);
    ceview.Lock();
    ce.Iterate();
    ceview.Unlock();
    long te = timer_us(timer);
    cout << "Iteration #" << i << " took: " << te << " us." << endl;
    cout << "Min Cost=" << ce.J << "\tAve Cost=" << ce.ce.Jave << endl;
    getchar();
  }

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
  viewer->frameName = "videos/sys";

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
