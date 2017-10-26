#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "gcarview.h"
#include "utils.h"
#include "se2.h"
#include "lqcost.h"
#include "pqpdem.h"
#include "constraintcost.h"
#include "multicost.h"
#include "demview.h"
#include "unistd.h"


using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<GcarState, 4, 2> GcarDdp;
typedef PqpDem<GcarState, 4, 2> GcarPqpDem;
typedef ConstraintCost<GcarState, 4, 2, Dynamic, 1> PqpDemCost;

int GcarToBody3dState(Body3dState &xb, const GcarState &xc)
{
  xb.p[0] = xc.g(0,2);
  xb.p[1] = xc.g(1,2);
  xb.p[2] = 0.5;
  return 2; 
}

bool DesiredState(GcarState &x, double t, double R, double v) 
{
  if (t < 0)
    return false;
  
  double w = v/R;   // angular velocity around circle
  double a = w*t;   // angle around track
  
  Vector3d q(a + M_PI/2, cos(a)*R, sin(a)*R);
  SE2::Instance().q2g(x.g, q); // set initial pose
  x.v = v;

  return true;
}

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(44.875, 46.875, -0.4, 0.050001, -20);

  double Rd = 5;   // desired circle radius
  double vd = 2;   // desired forward speed

  int N = 32;
  double tf = 3;
  double h = tf/N;

  SE2 &se2 = SE2::Instance();  
  Gcar sys;

  GcarState xf;

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<GcarState> xs(N+1);
  Vector3d q0(0.1, .1, .1);
  se2.q2g(xs[0].g, q0); // set initial pose
  xs[0].v = vd;

  // controls
  vector<Vector2d> us(N);
  for (int i = 0; i < N; ++i) {
     us[i].setZero();
     //   us[i] = Vector2d(0.01,.1);
     //       us[N/2+i] = Vector2d(0.01,-.1);
  }
  

  LqCost< GcarState, 4, 2> lqcost(sys, tf, xf);
  lqcost.Q.setZero();
  lqcost.R(0,0) = .005;
  lqcost.R(1,1) = .001;

  string mapName("maps/gcar_track1.ppm");
  double mapCellSize = 2*Rd/500;
  double mapHeightScale = 100;
  const double o[3] = {-Rd, -Rd, -mapHeightScale+2};
  //  params.GetString("mapName", mapName);
  //  params.GetDouble("mapCellSize", mapCellSize);
  //  params.GetDouble("mapHeightScale", mapHeightScale);
  Dem dem(mapName.c_str(), mapCellSize, mapHeightScale, o);
  //  dem.Convolve(1);

  DemView dv(dem);
  //  dv.wire = true;
  if (viewer)
    viewer->Add(dv);
  
  GcarPqpDem pqp(dem, .2, 0.5);
  pqp.func = GcarToBody3dState;
  PqpDemCost pqpcost(sys, tf, pqp);
  pqpcost.b = 1000;//2*pqpcost.b;

  MultiCost<GcarState, 4, 2> mcost(sys, tf);
  mcost.costs.push_back(&lqcost);

  //  Cylinder cyl(Vector3d(4,4,0), .5, 2, .1);
  //  CylinderView cview(cyl);
  
  //  viewer->Add(cview);
  mcost.costs.push_back(&pqpcost);

  GcarView view(sys, &xs, &us);
  viewer->Add(view);

  double t = h;

  vector<GcarState> xfs;
  DesiredState(xf, t, Rd, vd);
  xfs.push_back(xf);

  GcarView gview(sys, &xfs);
  gview.rgba[0] = 1;  gview.rgba[1] = 0;  gview.rgba[2] =0;
  viewer->Add(gview);

  getchar();

  struct timeval timer;


  while (t < 1000) {
    GcarDdp ddp(sys, mcost, ts, xs, us);
    // ddp.mu = 1;    
    // ddp.debug = false; // turn off debug for speed
      
    for (int i = 0; i < 20; ++i) {
      
      timer_start(timer);
      ddp.Iterate();
      long te = timer_us(timer);
      
      cout << "Iteration #" << i << " took: " << te << " us." << endl;    
    }
    getchar();  
    
    t = t + h;
    xs[0] = xs[1];
    rotate(us.begin(),us.begin()+1,us.end());
    us.back().setZero();
    //    for (int i = 0; i < N; ++i)
      //      us[i].setZero();

    // update desired state
    DesiredState(xf, t, Rd, vd);
    gview.Lock();
    xfs.push_back(xf);
    gview.Unlock();


  }
  //  for (int k = 0; k <= N; ++k)
  //    cout << ddp.xs[k] << "|" << endl;  
  //  cout << "xf=" << ddp.xs.back() << endl;
  
  
  cout << "done!" << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{

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
