#include <iomanip>
#include <iostream>
#include "sddp.h"
#include "viewer.h"
#include "rccarview.h"
#include "utils.h"
#include "rnlqcost.h"
#include "params.h"
#include "bulletrccar.h"
#include "bulletworld.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef SDdp<Vector4d, 4, 2> RccarDdp;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(-5, 51, -0.2, -0.15, -2.3);

  int N = 64;        // number of segments
  double tf = 5;    // time horizon

  int iters = 30;

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);

  params.GetInt("iters", iters);
  

  double h = tf/N;   // time step

  //Create Bullet world and rccar system:
  BulletWorld world(true);//Set the up axis as z for this world

  Bulletrccar sys(world);
  sys.kp_steer = 1.0f;

  //Load Ground
  {
    btCollisionShape *groundShape = world.CreateGroundPlane(50, 50);//20 by 20 long plane
    btTransform tr;
    tr.setOrigin(btVector3(0, 0, 0));
    tr.setRotation(btQuaternion(0,0,0));
    world.LocalCreateRigidBody(0,tr, groundShape);
  }

  //  sys.U.lb[1] = tan(-M_PI/5);
  //  sys.U.ub[1] = tan(M_PI/5);

  // initial state
  Vector4d x0(1,1,0,0);
  params.GetVector4d("x0", x0);

  // final state
  Vector4d xf(0,0,0,0);
  params.GetVector4d("xf", xf);  

  // cost
  RnLqCost<4, 2> cost(sys, tf, xf);
  VectorXd Q(4);
  if (params.GetVectorXd("Q", Q))
    cost.Q = Q.asDiagonal();
  
  VectorXd Qf(4);
  if (params.GetVectorXd("Qf", Qf))
    cost.Qf = Qf.asDiagonal();
  
  VectorXd R(2);
  if (params.GetVectorXd("R", R)) 
    cost.R = R.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Vector4d> xs(N+1);
  // initial state
  xs[0] = x0;

  // initial controls
  vector<Vector2d> us(N);
  for (int i = 0; i < N/2; ++i) {
    us[i] = Vector2d(.1, .05);
    us[N/2+i] = Vector2d(.1, -.05);
  }
  us[N-1] = Vector2d(0,0);
  
  RccarDdp ddp(sys, cost, ts, xs, us);  
  ddp.mu = .01;
  ddp.dxscale<<0.01,0.01,0.005,0.05;
  params.GetVector2d("scale_du",(ddp.duscale));
  params.GetDouble("mu", ddp.mu);

  RccarView view(sys, &ddp.xs);
  
  viewer->Add(view);

  struct timeval timer;
  // ddp.debug = false; // turn off debug for speed
  //getchar();

    timer_start(timer);
  for (int i = 0; i < iters; ++i) {
    ddp.Iterate();
    //getchar();
  }

  long te = timer_us(timer);
  cout << "Iterations" << iters << " took: " << te << " us." << endl;
  getchar();
  
  cout << xs[N] << endl;

  //  xs[1][3]  velocity
  //atan(us[0][1]) steering angle
 
  cout << "done!" << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{

  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("../../bin/rccar.cfg"); 


#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "../../logs/rccar/frames/frame";

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
