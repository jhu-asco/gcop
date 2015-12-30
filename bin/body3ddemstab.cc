#include <iomanip>
#include <iostream>
//#include "body3dcontroller.h"
#include "body3davoidcontroller.h"
//#include "gavoidcontroller.h"
#include "viewer.h"
#include "body3dview.h"
#include "utils.h"
#include "demview.h"
#include "pqpdem.h"
#include "params.h"
#include "dsl/gridsearch.h"
#include "dsl/gridcost.h"
#include "dsl/grid2d.h"
#include "dsl/grid2dconnectivity.h"

/*
  Tests Gyroscopic obstacle avoidance for a rigid body in a geometric terrain
 */


using namespace std;
using namespace Eigen;
using namespace gcop;
using namespace dsl;

typedef PqpDem<Body3dState, 12, 6> Body3dPqpDem;

Params params;


// produce a geometric trajectory b/n x0 and xf on a dem using grid D*-Lite
void GetTraj(vector<Body3dState> &gxs, Dem &dem, GridSearch<2> &gdsl,
   const Body3dState &x0, const Body3dState &xf) {
  
  int i0,j0,ig,jg;
  dem.Point2Index(i0, j0, x0.p[0], x0.p[1]);
  dem.Point2Index(ig, jg, xf.p[0], xf.p[1]);
  gdsl.SetStart(Vector2d(j0, i0));
  gdsl.SetGoal(Vector2d(jg, ig));
  GridPath<2> path, optPath;
  gdsl.Plan(path);
  gdsl.OptPath(path, optPath, 2);
  for (int i = 0; i < optPath.cells.size(); ++i) {
    Body3dState x = xf;
    dem.Index2Point(x.p[0], x.p[1], optPath.cells[i].c[1], optPath.cells[i].c[0]);
    x.p[2] = x0.p[2];
    
    // if not last point
    Vector3d v(0,0,0);
    if (i < optPath.cells.size() - 1) {
      Vector3d pa;
      Vector3d pb;
      dem.Index2Point(pa[0], pa[1], optPath.cells[i].c[1], optPath.cells[i].c[0]);
      dem.Index2Point(pb[0], pb[1], optPath.cells[i+1].c[1], optPath.cells[i+1].c[0]);
      pa[2] = x0.p[2];
      pb[2] = x0.p[2];
      Vector3d v = pb - pa;
      v = v/v.norm();
      v = v*20;
    }
    x.v = v;      
    gxs.push_back(x);
  }
}


void solver_process(Viewer* viewer)
{

 //#define MAP_SIMPLE  -- uncomment to use a simple terrain

  Body3d<> sys;
  int N = 500;      // discrete trajectory segments
  double tf = 10;   // time-horizon

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);

#ifdef MAP_SIMPLE

  if (viewer)
    viewer->SetCamera(4.5, 45, -10, -8, -10);
  Dem dem("maps/simple.ppm", .5, 15);
  dem.Convolve(2);

  // cost 
  Body3dState xf;
  xf.Clear();
  xf.p << 5, 5, 5;

  Body3dPqpDem pqp(dem, .1);  
  Body3dAvoidController<> ctrl(sys, &xf, 0, &pqp);
  ctrl.avoidCtrl->k = 5; // avoidance gain

  // states
  vector<Body3dState> xs(N+1);
  xs[0].Clear();
  xs[0].p << 20, 17, 10;
  
#else

  // cost 
  Body3dState x0;
  x0.Clear();
  x0.p << 46, 82, 2;

  Body3dState xf;
  xf.Clear();
  xf.p << 160, 125, 2;

  VectorXd qv0(12);
  params.GetVectorXd("x0", qv0);  
  SO3::Instance().q2g(x0.R, qv0.head(3));    
  x0.p = qv0.segment<3>(3); x0.w = qv0.segment<3>(6); x0.v = qv0.tail<3>(); 

  VectorXd qvf(12);
  params.GetVectorXd("xf", qvf);  
  SO3::Instance().q2g(xf.R, qvf.head(3));    
  xf.p = qvf.segment<3>(3); xf.w = qvf.segment<3>(6); xf.v = qvf.tail<3>(); 

  vector<Body3dState> xs(N+1);
  xs[0] = x0;
  
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

  params.GetVector6d("kp", ctrl.stabCtrl.Kp);  
  params.GetVector6d("kd", ctrl.stabCtrl.Kd);  
  params.GetDouble("ko", ctrl.avoidCtrl->k); 
  params.GetDouble("kb", ctrl.avoidCtrl->kb);  
  params.GetDouble("sr", ctrl.avoidCtrl->sr);

#endif

  DemView dv(dem);
  if (viewer)
    viewer->Add(dv);

  double h = tf/N; // time-step

  // trajectory view
  Body3dView<6> view(sys, &xs);
  view.rgba[0] = 0; view.rgba[1] = 0; view.rgba[2] = 1;

  // goal view
  vector<Body3dState> xfs;
  xfs.push_back(xf);
  Body3dView<6> goalView(sys, &xfs);
  goalView.rgba[0] = 1;

  if (viewer) {
    viewer->Add(view);  
    viewer->Add(goalView);  
  }
  double temp[50000];
  double temp2[50000];

  vector<Body3dState> gxs;
  Grid2d grid(dem.nj, dem.ni, dem.data, 1, 1, 1, 1);
  Grid2dConnectivity grid_connectivity(grid);
  GridCost<2> gridcost;
  GridSearch<2> gdsl(grid, grid_connectivity, gridcost);
  GetTraj(gxs, dem, gdsl, x0, xf);
  
  Body3dView<6> dslView(sys, &gxs);
  if (viewer)
    viewer->Add(dslView);  

  Vector6d u; 

  sys.U.bnd = true;
  sys.U.lb.setConstant(-50);
  sys.U.ub.setConstant(50);

  int j = 0;
  for (int i = 0; i < N; ++i) {

    ctrl.stabCtrl.xd = &gxs[j];    
    double t = i*h;
    ctrl.Set(u, t, xs[i]);
    sys.U.Clip(u);
    sys.Step(xs[i+1], t, xs[i], u, h);  
    
    // if close to current waypoint, then move to next
    Vector3d d = xs[i+1].p - ctrl.stabCtrl.xd->p;
    if (d.norm() < 20 && j < gxs.size()-1)
      ++j;
  }

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
  viewer->frameName = "body3d/frames/frame";

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
