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
#include "params.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Ddp<Body3dState, 12, 6> Body3dDdp;
typedef PqpDem<Body3dState, 12, 6> Body3dPqpDem;
typedef ConstraintCost<Body3dState, 12, 6, Dynamic, 1> PqpDemCost;

Params params;


void body3d2quad(Body3dState &x, Vector4d &u,
                 const Body3dState &x0, const Vector6d &u0,
                 const Body3d<> &sys,
                 bool gr = true)
{
  u.head(3) = u0.head(3);   
  if (gr)
    u[3] = (x0.first*u0.tail(3)).norm();
  else      
    u[3] = (x0.first*u0.tail(3) + Vector3d(0,0,9.81*sys.m)).norm();
  
  // required control force in spatial frame 
  Vector3d f = x0.first*u0.tail(3);
  if (!gr)
    f += Vector3d(0,0,9.81*sys.m);
  
  Vector3d bz = f/f.norm();
  Vector3d by = bz.cross(Vector3d(1,0,0));
  Vector3d bx = by.cross(bz);

  x = x0;
  x.first.col(0) = bx;
  x.first.col(1) = by;
  x.first.col(2) = bz;
}

void solver_process(Viewer* viewer)
{
  
  int N = 96;      // discrete trajectory segments
  double tf = 3;   // time-horizon
  int iters = 5; 

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);
  params.GetInt("iters", iters);

  double h = tf/N; // time-step
  
  double mass = 2;
  params.GetDouble("mass", mass);
  // system
  Body3d<> sys(Vector3d(.6, .6, .2), mass);
  bool gr = true;

  double cr = 0;
  params.GetDouble("cr", cr);

  if (gr)
    sys.fp << 0, 0, -9.81*sys.m;

  sys.U.bnd = true;
  params.GetBool("bnd", sys.U.bnd);
  sys.U.bndType = sys.U.BND_ELLIPSOID;
  double umax = 30;
  params.GetDouble("umax", umax);
  double wi = 1/(umax*umax);
  sys.U.w << 0, 0, 0, wi, wi, wi; // ellispoidal bound weights

  // states
  vector<Body3dState> xs(N+1);

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

  DemView dv(dem);
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
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // initial controls (e.g. hover at one place)
  vector<Vector6d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].setZero();
    us[i].tail<3>() = -sys.fp;
  }

  Body3dView<6> view(sys, &xs);
  view.rgba[0] = 1; view.rgba[1] = 0; view.rgba[2] = 0;
  if (viewer)
    viewer->Add(view);  

  double temp[50000];

  Body3dPqpDem pqp(dem, cr);

  double temp2[50000];

  PqpDemCost pqpcost(sys, tf, pqp);
  
  MultiCost<Body3dState, 12, 6> mcost(sys, tf);
  mcost.costs.push_back(&cost);   
  mcost.costs.push_back(&pqpcost);
    
  Body3dDdp ddp(sys, mcost, ts, xs, us);
  ddp.eps = 1e-3;
  ddp.mu = 10;
  params.GetDouble("mu", ddp.mu);
  params.GetDouble("eps", ddp.eps);


  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed

  //getchar();

  vector<Body3dState> hxs;
  vector<Vector4d> hus;
  vector<double> hts;
  hxs.push_back(xs[0]);
  hts.push_back(0);
  Qrotor hrotor;
  QrotorView hv(hrotor, &hxs);
  hv.lineWidth = 5;
  hv.rgba[0] = 0;   hv.rgba[1] = 0;   hv.rgba[2] = 1;
  if (viewer)
    viewer->Add(hv);

  bool xfStop = true; // stop when predicted rather than current is near xf

  while(1) {    
    struct timeval timer;
    //  ddp.debug = false; // turn off debug for speed
    
    //    pqpcost.b = 1;
    timer_start(timer);
    for (int i = 0; i < iters; ++i) {
      ddp.Iterate();
      pqpcost.b = 1000;//2*pqpcost.b;
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
    Vector4d hu;
    hus.push_back(hu);

    body3d2quad(hxs.back(), hus.back(), xs[0], us[0], sys);

    if (xfStop) {
      if ((xs.back().second.head(3) - xf.second.head(3)).norm() < 5) {
      for (int k=0; k <xs.size(); ++k) {
        hxs.push_back(xs[k]);
          hts.push_back(hts.back()+h);
          Vector4d hu;
          hus.push_back(hu);          
          body3d2quad(hxs.back(), hus.back(), xs[k], us[k], sys);          
        }
        
        break;
      }
    } else {
      if ((xs[0].second.head(3)-xf.second.head(3)).norm() < 5)
        break;
    }
  }

  fstream fstr;
  fstr.open("traj.txt", std::ios::out | std::ios::trunc);
  fstr.precision(20);
  //  fstr << hxs.size()-1 << " 0" << endl;
  for (int i = 0; i < hxs.size(); ++i) {
    const Vector3d &p = hxs[i].second.head<3>();
    fstr << hts[i] << " " << hxs[i].second.tail<6>().transpose() << " ";

    Vector3d rpy;
    SO3::Instance().g2q(rpy, hxs[i].first);

    fstr << rpy[0] << " " << rpy[1] << " " << rpy[2] << " ";
    
    fstr << p[0] << " " << p[1] << " " << p[2] << " ";
    if (i < hus.size())
      fstr << hus[i].transpose() << endl;
    else 
      fstr << "0 0 0 0" << endl;
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
  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("../../bin/body3drhc.cfg");

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

