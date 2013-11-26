#include <iomanip>
#include <iostream>
#include "dmoc.h"
#include "viewer.h"
#include "airplaneview.h"
#include "airplanecost.h"
#include "utils.h"
#include "se2.h"
#include "normal2dview.h"
#include <cmath>

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<pair<Matrix3d, Vector2d>, 5, 2> GunicycleDmoc;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(25, 27, -3.3, -1.1, -3.65);

  int N = 128;
  double tf = 10;
  double h = tf/N;

  SE2 &se2 = SE2::Instance();  
  Gunicycle sys;

  M3V2d xf(se2.Id, Vector2d(0, 1));
  Vector3d qf(0, 10, 0);
  se2.q2g(xf.first, qf);

  AirplaneCost cost(tf, xf, N);  

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<pair<Matrix3d, Vector2d> > xs(N+1);
  Vector3d q0(0, 0, 0);
  se2.q2g(xs[0].first, q0);
  xs[0].second[0] = 0;
  xs[0].second[1] = 1;

  // controls
  vector<Vector2d> us(N);
  for (int i = 0; i < N/2; ++i) {
    //    us[i] = Vector2d(0.01,.1);
    //    us[N/2+i] = Vector2d(0.01,-.1);
  }  

  GunicycleDmoc dmoc(sys, cost, ts, xs, us);
  dmoc.mu = 0;

  AirplaneView view(sys, &dmoc.xs);
  viewer->Add(view);

  //  double cs = .1;
  //  Dem dem(10, 10, cs, 1);
  //  Vector2d po(5, 5);
    //  dem.Set(10, 10, 5);

  //  DemView dv(dem);
  Normal2dView dv(&cost.gds);
  dv.wire = true;
  dv.s = 20;
  viewer->Add(dv);

  struct timeval timer;

  dmoc.debug = false; // turn off debug for speed

  //  for (cost.n = 0; cost.n < 5; cost.n += .5) {
  for (double s = 2; s >= .01; s -= .1) {
    getchar();
    
    cost.SetObstacle(s, 125);
    
    for (int i = 0; i < 10; ++i) {      
      timer_start(timer);
      dmoc.Iterate();
      long te = timer_us(timer);      
      cout << "Iteration #" << i << " took: " << te << " us." << endl;    
    }
  }
  //  for (int k = 0; k <= N; ++k)
  //    cout << dmoc.xs[k] << "|" << endl;  
  //  cout << "xf=" << dmoc.xs.back() << endl;  
  
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
  viewer->frameName = "../../logs/airplane/frames/frame";

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
