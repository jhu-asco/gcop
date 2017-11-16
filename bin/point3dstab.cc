#include <iostream>
#include "viewer.h"
#include "point3dview.h"
#include "point3dcontroller.h"
#include "utils.h"
#include "disk.h"
#include <unistd.h>

using namespace std;
using namespace Eigen;
using namespace gcop;

//typedef Ddp<Point3dState, 4, 2> Point3dDdp;
//typedef Disk<Point3dState, 4, 2> Point3dDisk;
//typedef ConstraintCost<Point3dState, 4, 2, Dynamic, 1> DiskCost;


void solver_process(Viewer* viewer)
{
  viewer->SetCamera(22.625, 45, 1.25, 3.15, -10);

  int N = 500;
  double tf = 10;
  double h = tf/N;

  Point3d sys;

  Point3dState xf;

  // states
  vector<Point3dState> xs(N+1);
  xs[0].q << -5, -5, 2;
  xs[0].v << 0, 1, 0;
  
  //  Point3dDisk disk(Vector3d(-2.5,-2.5), 2, 0);
  //  DiskCost dcost(sys, tf, disk);
  
  Point3dController ctrl(sys);

  Point3dView view(sys, &xs);
  viewer->Add(view);  

  Vector3d u; 

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
