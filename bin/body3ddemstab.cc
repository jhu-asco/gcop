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

/*
  Tests Gyroscopic obstacle avoidance for a rigid body in a geometric terrain
 */


using namespace std;
using namespace Eigen;
using namespace gcop;

typedef PqpDem<Body3dState, 12, 6> Body3dPqpDem;

Params params;

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
  Body3dState xf(Matrix3d::Identity(), Vector9d::Zero());
  xf.second[0] = 5;  
  xf.second[1] = 5;  
  xf.second[2] = 5;

  Body3dPqpDem pqp(dem, .1);  
  Body3dAvoidController<> ctrl(sys, &xf, 0, &pqp);
  ctrl.avoidCtrl->k = 5; // avoidance gain

  // states
  vector<Body3dState> xs(N+1);
  xs[0].first.setIdentity();
  xs[0].second.setZero();
  xs[0].second[0] = 20;
  xs[0].second[1] = 17;
  xs[0].second[2] = 10;
  
#else

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

  Vector6d u; 
  for (int i = 0; i < N; ++i) {
    double t = i*h;
    ctrl.Set(u, t, xs[i]);
    sys.Step(xs[i+1], t, xs[i], u, h);    
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
