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

#include "hrotorview.h"


using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<Body3dState, 12, 6> Body3dDdp;
typedef PqpDem<Body3dState, 12, 6> Body3dPqpDem;
typedef ConstraintCost<Body3dState, 12, 6, Dynamic, 1> PqpDemCost;


void solver_process(Viewer* viewer)
{

  int N = 96;      // discrete trajectory segments
  double tf = 3;   // time-horizon
  double h = tf/N; // time-step

  // system
  Body3d<> sys;

  // states
  vector<Body3dState> xs(N+1);

  // cost 
  Body3dState xf(Matrix3d::Identity(), Vector9d::Zero());


#define MAP_CITY

#ifdef MAP_CITY
  if (viewer)
    viewer->SetCamera(90, 34, 100.2, -78, -32);
  
  Dem dem("maps/pic2_dilated.ppm", .25, 30);
  //  dem.Dilate(5);
  //  dem.Convolve(.3);

    xs[0].first.setIdentity();
  xs[0].second[0] = 46;
  xs[0].second[1] = 82;
  xs[0].second[2] = 2;

  xf.second[0] = 165;  
  xf.second[1] = 130;  
  xf.second[2] = 2;

#else
  if (viewer)
    viewer->SetCamera(4.5, 45, -10, -8, -10);
  Dem dem("maps/simple.ppm", .5, 15);
  dem.Convolve(2);

  xs[0].first.setIdentity();
  xs[0].second[0] = 20;
  xs[0].second[1] = 17;
  xs[0].second[2] = 10;

  xf.second[0] = 5;  
  xf.second[1] = 5;  
  xf.second[2] = 3;

#endif

  DemView dv(dem);
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
  vector<Vector6d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].setZero();
  }

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
    
  Body3dDdp ddp(sys, mcost, ts, xs, us);
  ddp.eps = 1e-3;
  ddp.mu = 10;

  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed

  getchar();

  vector<Body3dState> hxs;
  vector<double> hts;
  hxs.push_back(xs[0]);
  hts.push_back(0);
  Hrotor hrotor;
  HrotorView hv(hrotor, &hxs);
  if (viewer)
    viewer->Add(hv);

  while(1) {    
    struct timeval timer;
    //  ddp.debug = false; // turn off debug for speed
    
    //    pqpcost.b = 1;
    timer_start(timer);
    for (int i = 0; i < 4; ++i) {
      ddp.Iterate();
      pqpcost.b = 5000;//2*pqpcost.b;
    }
    long te = timer_us(timer);
    cout << "Took " << te << " us." << endl;        
    //    getchar();
    xs[0] = xs[1];
    std::rotate(us.begin(), us.begin() + 1, us.end());
    us.back() = us[us.size()-2];
    
    for (int i = 0; i < N; ++i) {
      //      us[i].setZero();
    }

    ddp.Update();
    cout << xs[0].first << " " << xs[0].second  << endl;
    cout << "Moved forward" << endl;
    //    getchar();
    hxs.push_back(xs[0]);
    hts.push_back(hts.back()+h);

    // required control force in spatial frame 
    Vector3d f = xs[0].first*us[0].tail(3) + Vector3d(0,0,9.81*sys.m);
    Vector3d bz = f/f.norm();
    Vector3d by = bz.cross(Vector3d(1,0,0));
    Vector3d bx = by.cross(bz);
    hxs.back().first.col(0) = bx;
    hxs.back().first.col(1) = by;
    hxs.back().first.col(2) = bz;

    if ((xs[0].second.head(3)-xf.second.head(3)).norm() < 2)
      break;
  }

  fstream fstr;
  fstr.open("traj.txt", std::ios::out | std::ios::trunc);
  fstr.precision(20);
  fstr << hxs.size()-1 << " 0" << endl;
  for (int i = 0; i < hxs.size(); ++i) {
    const Vector3d &p = hxs[i].second.head<3>();
    fstr << hts[i] << " 0 0 0 0 0 0 ";

    Vector3d rpy;
    SO3::Instance().g2q(rpy, hxs[i].first);

    fstr << rpy[0] << " " << rpy[1] << " " << rpy[2] << " ";
    fstr << p[0] << " " << p[1] << " " << p[2] << " 0 0 0 0" << endl;
  }
  fstr.close();

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

