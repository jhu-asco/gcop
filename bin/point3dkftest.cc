#include <iostream>
#include "utils.h"
#include "so3.h"
#include "kalmanpredictor.h"
#include "kalmancorrector.h"
#include "unscentedpredictor.h"
#include "unscentedcorrector.h"
#include "point3d.h"
#include "point3dgps.h"

using namespace std;
using namespace gcop;
using namespace Eigen;

typedef KalmanPredictor<Point3dState, 6, 3, Dynamic> Point3dKalmanPredictor;
typedef KalmanCorrector<Point3dState, 6, 3, Dynamic, Vector3d, 3> Point3dGpsKalmanCorrector;

int main(int argc, char** argv)
{
  Point3d point3d;

  Point3dGps<> gps;

  Point3dKalmanPredictor kp(point3d);
  Point3dGpsKalmanCorrector kc(point3d.X, gps);

  int N = 1000;

  vector<Point3dState> xts(N);   // true states
  xts[0].v << .1, 0, .2;         // true (initially unknown) velocity

  vector<Point3dState> xs(N);   // estimated trajectory
  xs[0].P.topLeftCorner<3,3>().diagonal().setConstant(.1);  // q
  xs[0].P.bottomRightCorner<3,3>().diagonal().setConstant(.1);  // v

  double dt = .01;

  struct timeval timer;
  
  for (int i = 0; i < N-1; ++i) {
    double t = i*dt;

    // Constant velocity model no acceleration
    Vector3d u(0, 0, 0); // no inputs

    // generate true
    point3d.Step(xts[i+1], t, xts[i], u, dt);

    Point3dState x;

    timer_start(timer);
    kp.Predict(x, t, xs[i], u, dt);
    long us = timer_us(timer);
    cout << "Predict took " << us << " us." << endl;
    cout << "qt: " << xts[i].q.transpose() << endl;
    cout << "qe: " << x.q.transpose() << endl;

    // noisy measurements of position
    Vector3d z = xts[i+1].q + Vector3d(gps.sxy*randn(), gps.sxy*randn(), gps.sz*randn());
    cout << "z=" << z.transpose() << endl;

    timer_start(timer);
    kc.Correct(xs[i+1], t, x, u, z);
    us = timer_us(timer);
    cout << "Update took " << us << " us." << endl;

    cout << "True position: " << xts[i+1].q << endl;
    cout << "Estim position: " << xs[i+1].q << endl;

    cout << "True velocity: " << xts[i+1].v << endl;
    cout << "Estim velocity: " << xs[i+1].v << endl;

    cout << "P=: " << xs[i+1].P << endl;    
  }
}
