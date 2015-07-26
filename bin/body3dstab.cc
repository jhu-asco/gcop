#include <iostream>
#include "viewer.h"
#include "body3dview.h"
#include "body3dcontroller.h"
#include "utils.h"
#include "disk.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

//typedef Ddp<Body3dState, 4, 2> Body3dDdp;
//typedef Disk<Body3dState, 4, 2> Body3dDisk;
//typedef ConstraintCost<Body3dState, 4, 2, Dynamic, 1> DiskCost;



void solver_process(Viewer* viewer)
{
  viewer->SetCamera(-7.25, 61, -2.8, -2.55, -10);

  int N = 500;
  double tf = 10;
  double h = tf/N;

  Body3d<6> sys;

  Body3dState xf;
  xf.Clear();

  // states
  vector<Body3dState> xs(N+1);
  Vector3d e0(1.2, -2, 1);
  xs[0].Clear();
  SO3::Instance().exp(xs[0].R, e0);
  xs[0].p << 5, 5, 5;
  
  //  Body3dDisk disk(Vector3d(-2.5,-2.5), 2, 0);
  //  DiskCost dcost(sys, tf, disk);
  
  Body3dController<> ctrl(sys);

  Body3dView<> view(sys, &xs);
  viewer->Add(view);  

  Vector6d u; 

  for (int i = 0; i < N; ++i) {
    double t = i*h;
    ctrl.Set(u, t, xs[i]);
    sys.Step(xs[i+1], t, xs[i], u, h);
  }

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
