#include <iostream>
#include "utils.h"
#include "so3.h"
#include "kalmanpredictor.h"
#include "kalmancorrector.h"
#include "body3dsensor.h"

using namespace std;
using namespace gcop;
using namespace Eigen;

typedef KalmanPredictor<Body3dState, 12, 6, Dynamic> Pose3dKalmanPredictor;
typedef KalmanCorrector<Body3dState, 12, 6, Dynamic, Pose3d, 6> Pose3dKalmanCorrector;

int main(int argc, char** argv)
{
  Body3d<6> sys;
  sys.constantVelocity = true;
  Body3dSensor<6> sensor;

  Pose3dKalmanPredictor kp(sys);
  Pose3dKalmanCorrector kc(sys.X, sensor);

  int N = 1000;

  vector<Body3dState> xts(N);   // true states
  xts[0].Clear();
  xts[0].v << .1, 0, .2;         // true (initially unknown) velocity
  xts[0].w << .1, .05, .2;         // true (initially unknown) velocity

  vector<Body3dState> xs(N);   // estimated trajectory
  xs[0].Clear();

  xs[0].P.topLeftCorner<3,3>().diagonal().setConstant(.1);  // R
  xs[0].P.block<3,3>(3,3).diagonal().setConstant(.1);  // p
  xs[0].P.block<3,3>(6,6).diagonal().setConstant(.1);  // w
  xs[0].P.bottomRightCorner<3,3>().diagonal().setConstant(.1);  // v


  double dt = .01;

  struct timeval timer;
  
  for (int i = 0; i < N-1; ++i) {
    double t = i*dt;

    // Constant velocity model no acceleration
    Vector6d u; u.setZero();

    // generate true
    sys.Step(xts[i+1], t, xts[i], u, dt);

    Body3dState x;

    timer_start(timer);
    kp.Predict(x, t, xs[i], u, dt);
    long us = timer_us(timer);
    cout << "Predict took " << us << " us." << endl;
    cout << "qt: " << xts[i].p.transpose() << endl;
    cout << "qe: " << x.p.transpose() << endl;

    // noisy measurements of position
    Pose3d z;
    z.p = xts[i+1].p + sensor.sp.cwiseProduct(Vector3d(randn(), randn(), randn()));

    Matrix3d dR;
    SO3::Instance().exp(dR, sensor.sR.cwiseProduct(Vector3d(randn(), randn(), randn())));
    z.R = xts[i+1].R*dR;

    timer_start(timer);
    kc.Correct(xs[i+1], t, x, u, z);
    //xs[i+1] = x;
    us = timer_us(timer);
    cout << "Update took " << us << " us." << endl;

    cout << "True position: " << xts[i+1].p << endl;
    cout << "Estim position: " << xs[i+1].p << endl;
    
    Vector3d rpyt, rpy;
    SO3::Instance().g2q(rpyt, xts[i+1].R);
    SO3::Instance().g2q(rpy, xs[i+1].R);

    cout << "True rot: " << rpyt << endl;
    cout << "Estim rot: " << rpy<< endl;

    cout << "True w: " << xts[i+1].w << endl;
    cout << "Estim w: " << xs[i+1].w << endl;

    cout << "True velocity: " << xts[i+1].v << endl;
    cout << "Estim velocity: " << xs[i+1].v << endl;


    cout << "P=: " << xs[i+1].P << endl;    
  }
}
