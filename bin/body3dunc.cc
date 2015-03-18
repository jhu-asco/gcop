#include <iomanip>
#include <fstream>
#include "ddp.h"
#include "viewer.h"
#include "body3dview.h"
#include "utils.h"
#include "body3dcost.h"
#include "demview.h"
#include "pqpdem.h"
#include "constraintcost.h"
#include "multicost.h"

#include "qrotorview.h"


using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<Body3dState, 12, 6> Body3dDdp;
typedef PqpDem<Body3dState, 12, 6> Body3dPqpDem;
typedef ConstraintCost<Body3dState, 12, 6, Dynamic, 1> PqpDemCost;


void solver_process(Viewer* viewer)
{

  // phi=78.25 theta=55.875 x=66.9499 y=-69.1997 z=78.3501
  int N = 128;      // discrete trajectory segments
  double tf = 12;   // time-horizon
  double h = tf/N; // time-step

  // system
  Body3d<> sys;

  // states
  vector<Body3dState> xs(N+1);

  // cost 
  Body3dState xf(Matrix3d::Identity(), Vector9d::Zero());


  if (viewer)
    viewer->SetCamera(90, 34, 100.2, -78, -32);
  
  Dem dem("maps/pic2_dilated.ppm", .25, 30);
  Dem dem2("maps/pic2.ppm", .25, 30);  // for display purposes only

  //  dem.Dilate(5);
  //  dem.Convolve(.3);
  Body3dState x0;

  x0.first.setIdentity();
  x0.second[0] = 95;
  x0.second[1] = 87;
  x0.second[2] = 3;

  //  xf.second[0] = 165;  
  //  xf.second[1] = 130;  
  xf.second[0] = 110;
  xf.second[1] = 120;  

  xf.second[2] = 3;

  xs[0] = x0;


  DemView dv(dem2);
  //  dv.wire = true;
  if (viewer)
    viewer->Add(dv);


  Body3dCost<6> cost(sys, tf, xf);  

  double s = .01;
  cost.Qf(0,0) = 2; cost.Qf(1,1) = 2; cost.Qf(2,2) = 2;
  cost.Qf(3,3) = 20; cost.Qf(4,4) = 20; cost.Qf(5,5) = 20;

  //  cost.Q(0,0) = 1; cost.Q(1,1) = 1; cost.Q(2,2) = 1;
  //  cost.Q(3,3) = 2; cost.Q(4,4) = 2; cost.Q(5,5) = 2;
  cost.Q(11,11) = 10;

  cost.Qf(6,6) = 5; cost.Qf(7,7) = 5; cost.Qf(8,8) = 5;
  cost.Qf(9,9) = 50; cost.Qf(10,10) = 50; cost.Qf(11,11) = 50;

  cost.Q = s*cost.Q;
  cost.Qf = s*cost.Qf;
  
  cost.R(0,0) = .01; cost.R(1,1) = .01; cost.R(2,2) = .01; 
  cost.R(3,3) = .1; cost.R(4,4) = .1; cost.R(5,5) = .1;

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // initial controls (e.g. hover at one place)

  Body3dView<6> view(sys, &xs);
  if (viewer)
    viewer->Add(view);  

  double temp[50000];

  Body3dPqpDem pqp(dem, .01);

  double temp2[50000];

  PqpDemCost pqpcost(sys, tf, pqp);
  
  MultiCost<Body3dState, 12, 6> mcost(sys, tf);
  mcost.costs.push_back(&cost);   
  mcost.costs.push_back(&pqpcost);
    
  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed

  getchar();

  vector<Vector4d> hus(N);

  Qrotor hrotor;

  int M = 20;

  vector<Body3dState> trajs[M];
  QrotorView *views[M];

  VectorXd Vs(M);
  
  for (int j = 0; j <M; ++j) {    
    sys.fp << .2*(RND-.5), .2*(RND-.5), 0;
    
    xs[0].second.head(3) = x0.second.head(3) + 2*Vector3d(RND, RND, RND);

    vector<Vector6d> us(N);
    for (int i = 0; i < N; ++i) {
      us[i].setZero();
    }
    
    Body3dDdp ddp(sys, mcost, ts, xs, us);
    ddp.eps = 1e-3;
    ddp.mu = 1;    

    trajs[j].resize(N);
    views[j] = new QrotorView(hrotor, &trajs[j]);
    views[j]->rgba[0] = 1;
    views[j]->rgba[1] = 0;
    views[j]->rgba[2] = 0;

    timer_start(timer);
    pqpcost.b = 1000;
    for (int i = 0; i < 20; ++i) {
      ddp.Iterate();
      //      pqpcost.b = 2*pqpcost.b;
    }
    long te = timer_us(timer);
    cout << "j=" << j << "  took " << te << " us." << endl;                
    Vs(j) = ddp.V;
    
    for (int  i = 0; i < N; ++i) {
      Vector4d hu;
      hu.head(3) = us[0].head(3);   
      hu[3] = (xs[0].first*us[0].tail(3) + Vector3d(0,0,9.81*sys.m)).norm();
      //  hus[i] = hu;
      
      // required control force in spatial frame 
      Vector3d f = xs[0].first*us[0].tail(3) + Vector3d(0,0,9.81*sys.m);
      Vector3d bz = f/f.norm();
      Vector3d by = bz.cross(Vector3d(1,0,0));
      Vector3d bx = by.cross(bz);
      trajs[j][i] = xs[i];  
      trajs[j][i].first.col(0) = bx;
      trajs[j][i].first.col(1) = by;
      trajs[j][i].first.col(2) = bz;
    }
    
    //    
    if (viewer)
      viewer->Add(*views[j]);
    
  }

  cout <<Vs.transpose() << endl;

  viewer->Remove(view);

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
  viewer->frameName = "logs/body3d/frames/frame";
  viewer->displayName = "logs/body3d/display/frame";

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

