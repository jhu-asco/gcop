#include <iostream>
#include "utils.h"
#include "so3.h"
#include "kalmanfilter.h"
#include "unscentedfilter.h"
#include "imu.h"
#include "imusensor.h"

using namespace std;
using namespace gcop;
using namespace Eigen;

typedef KalmanFilter<ImuState, 9, 3, Dynamic, Matrix<double, 3, 1>, 3> ImuKF;
typedef UnscentedFilter<ImuState, 9, 3, Dynamic, Matrix<double, 3, 1>, 3> ImuUKF;

int main(int argc, char** argv)
{
  Imu imu;

  bool mag = false;
  //  int nz = mag ? 6 : 3;

  ImuSensor<3> sensor;
  //  sensor.sa = 0;  // no acceleration bias drift
  //  sensor.sv = 0;  // no acceleration bias drift

  ImuKF kf(imu, sensor);
  ImuUKF ukf(imu, sensor);
  
  int N = 150;

  vector<ImuState> xts(N);   // true states
  xts[0].bg << 0.05, 0, 0;   // bias

  vector<ImuState> xs(N);   // estimated trajectory
  xs[0].P.topLeftCorner<3,3>().diagonal().setConstant(.1);
  xs[0].P.block<3,3>(3,3).diagonal().setConstant(1e-2);
  xs[0].P.bottomRightCorner<3,3>().diagonal().setConstant(1e-10);

  // true rotation

  double dt = .01;

  struct timeval timer;
  
  Vector3d qt;
  Vector3d q;

  for (int i = 0; i < N-1; ++i) {
    double t = i*dt;

    // true input (drifted angular velocity)
    Vector3d wt(.2, .1, 0); // this is the undrifted
    wt = wt + xts[i].bg;    // drifted angular velocity

    // generate true
    imu.Step(xts[i+1], t, xts[i], wt, dt);

    // drifted noisy gyro reading
    Vector3d w = wt + imu.sv*Vector3d::Ones()*randn();

    ImuState x;

    timer_start(timer);
    kf.Predict(x, t, xs[i], w, dt);
    long us = timer_us(timer);
    cout << "Predict took " << us << " us." << endl;

    SO3::Instance().g2q(qt, xts[i+1].R);
    SO3::Instance().g2q(q, x.R);
    cout << "rpyt: " << qt.transpose() << endl;
    cout << "rpye: " << q.transpose() << endl;

    // true measurements of acc and mag
    Vector3d at = xts[i+1].R.transpose()*sensor.a0;
    Vector3d mt = xts[i+1].R.transpose()*sensor.m0;

    // noisy measurements
    /*
    Vector6d z;
    z.head<3>() = at + xts[i+1].second.tail<3>() + sensor.sra*Vector3d::Ones()*randn();
    z.head<3>().normalize();
    z.tail<3>() = mt + sensor.srm*Vector3d::Ones()*randn();
    z.tail<3>().normalize();
    */
    Vector3d z;
    z = at + xts[i+1].ba + sensor.sra*Vector3d::Ones()*randn();
    z.normalize();

    cout << "z=" << z.transpose() << endl;

    cout << "P-=" << x.P << endl;

    //    if (i==N-2) {
      timer_start(timer);
      kf.Update(xs[i+1], t, x, w, z);
      //xs[i+1] = x;
      us = timer_us(timer);
      cout << "Update took " << us << " us." << endl;
      //    }

    cout << "True attitude: " << xts[i+1].R << endl;
    cout << "Estim attitude: " << xs[i+1].R << endl;

    cout << "True biases: " << xts[i+1].bg.transpose() << endl;
    cout << "Estim biases: " << xs[i+1].bg.transpose() << endl;

    cout << "P=: " << xs[i+1].P << endl;
    
    if (std::isnan(xs[i+1].R(0,0)))
      break;
    
    SO3::Instance().g2q(qt, xts[i+1].R);
    SO3::Instance().g2q(q, xs[i+1].R);
    cout << "rpyt: " << qt.transpose() << endl;
    cout << "rpye: " << q.transpose() << endl;
    
    Vector3d dR;
    SO3::Instance().log(dR, xs[i+1].R.transpose()*xts[i+1].R);
    cout << "rotation error: " << dR.norm() << endl;
    //    cout << "L_2 norm (state) =" << (xt - ukf.x).norm() << endl;
    //    cout << "L_2 norm (quat) =" << (xt.head(4) - ukf.x.head(4)).norm() << endl;
  }
}
