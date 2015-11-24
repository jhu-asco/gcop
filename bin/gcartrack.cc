#include <iomanip>
#include <iostream>
#include "ddp.h"
#include "viewer.h"
#include "gcarview.h"
#include "utils.h"
#include "se2.h"
#include "lqcost.h"
#include "constraintcost.h"
#include "multicost.h"
#include "diskconstraint.h"
#include "diskview.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<GcarState, 4, 2> GcarDdp;
typedef DiskConstraint<GcarState, 4, 2> GcarDiskConstraint;
typedef ConstraintCost<GcarState, 4, 2, Dynamic, 1> DiskConstraintCost;


int GcarStateToVector2d(Vector2d &p, const GcarState &x)
{
  p[0] = x.g(0,2);   // x.g is the car SE(2) pose matrix 
  p[1] = x.g(1,2);
  return 1;  // index where position coordinates start (the car model coordinates are (theta,x,y,v)
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



bool MakeObstacles(vector<Disk> &disks, int N, double R, double dr) 
{
  if (R <= 0 || dr <=0)
    return false;
  for (int i = 0; i <N; ++i) {
    double a = RND*2*M_PI;
    disks.push_back(Disk(Vector2d(cos(a)*R, sin(a)*R), dr));
  }
  return true;
}


void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(44.875, 46.875, -0.4, 0.050001, -16);

  double Rd = 5;   // desired circle radius
  double vd = 2;   // desired forward speed

  int N = 32;
  double tf = 3;
  double h = tf/N;

  SE2 &se2 = SE2::Instance();  
  Gcar sys;

  GcarState xf; // final state

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

  MultiCost<GcarState, 4, 2> mcost(sys, tf);
  mcost.costs.push_back(&lqcost);
  
  // make obstacles and add as constarints
  vector<Disk> disks;
  MakeObstacles(disks, 6, Rd, .3);
  
  GcarDiskConstraint* constraints[disks.size()];
  DiskConstraintCost* dcosts[disks.size()];
  DiskView* dviews[disks.size()];

  for (int i = 0; i < disks.size(); ++i) {
    constraints[i] = new GcarDiskConstraint(disks[i], .3);
    constraints[i]->func = GcarStateToVector2d;    
    dcosts[i] = new DiskConstraintCost(sys, tf, *constraints[i]);
    dcosts[i]->b = 1000;
    mcost.costs.push_back(dcosts[i]);

    dviews[i] = new DiskView(disks[i]);
    viewer->Add(*dviews[i]);    
  }

  GcarView view(sys, &xs, &us);
  viewer->Add(view);

  double t = h;

  vector<GcarState> xfs;
  DesiredState(xf, t, Rd, vd);
  xfs.push_back(xf);

  GcarView gview(sys, &xfs);
  gview.rgba[0] = 1;  gview.rgba[1] = 0;  gview.rgba[2] =0;
  viewer->Add(gview);
  gview.renderSystem = false;

  // current (i.e. the system at the current time) view
  GcarView cview(sys);
  cview.x = &xs[0];
  cview.rgba[0] = 0;  cview.rgba[1] = 0;  cview.rgba[2] =1;
  viewer->Add(cview);

  getchar();

  struct timeval timer;

  while (t < 1000) {
    GcarDdp ddp(sys, mcost, ts, xs, us);
    // ddp.mu = 1;    
    ddp.debug = false; // turn off debug for speed

    for (int i = 0; i < 50; ++i) {      
      timer_start(timer);
      ddp.Iterate();
      long te = timer_us(timer);      
      //      cout << "Iteration #" << i << " took: " << te << " us." << endl;    
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
    
    // inefficient way to cut down the reference traj
    if (xfs.size() > xs.size())
      xfs.erase(xfs.begin());

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
