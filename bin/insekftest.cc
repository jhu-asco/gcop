#include <iostream>
#include "utils.h"
#include "so3.h"
#include "kalmanpredictor.h"
#include "kalmancorrector.h"
#include "unscentedpredictor.h"
#include "unscentedcorrector.h"
#include "ins.h"
#include "insimu.h"
#include "insgps.h"

using namespace std;
using namespace gcop;
using namespace Eigen;

typedef KalmanPredictor<InsState, 15, 6, Dynamic> InsKalmanPredictor;
typedef KalmanCorrector<InsState, 15, 6, Dynamic, Vector3d, 3> InsImuKalmanCorrector;
typedef KalmanCorrector<InsState, 15, 6, Dynamic, Vector3d, 3> InsGpsKalmanCorrector;


//typedef UnscentedPredictor<InsState, 9, 3, Dynamic> InsUnscentedPredictor;
//typedef UnscentedCorrector<InsState, 9, 3, Dynamic, Matrix<double, 3, 1>, 3> InsUnscentedCorrector;


int main(int argc, char** argv)
{
  Ins ins;

  // don't use magnetometer for now
  bool mag = false;
  //  int nz = mag ? 6 : 3;

  InsImu<3> imu;
  InsGps gps;

  InsKalmanPredictor kp(ins);
  InsImuKalmanCorrector kci(ins, imu);
  InsGpsKalmanCorrector kcg(ins, gps);

  //  InsUnscentedPredictor kp(ins);
  // InsUnscentedCorrector kc(ins, sensor);
  
  int N = 1000;

  vector<InsState> xts(N);   // true states
  xts[0].bg << 0.1, 0, 0;    // gyro bias
  //  xts[0].ba << 1e-5, 0, 0;    // acc bias
  xts[0].v << 1, 0, 0;       // true velocity

  vector<InsState> xs(N);   // estimated trajectory
  xs[0].P.topLeftCorner<3,3>().diagonal().setConstant(.1);  // R
  xs[0].P.block<3,3>(3,3).diagonal().setConstant(1e-2);     // bg
  xs[0].P.block<3,3>(6,6).diagonal().setConstant(1e-10);    // ba
  xs[0].P.block<3,3>(9,9).diagonal().setConstant(.01);      // p
  xs[0].P.block<3,3>(12,12).diagonal().setConstant(.04);    // v

  double dt = .01;

  struct timeval timer;
  
  Vector3d qt;
  Vector3d q;

  for (int i = 0; i < N-1; ++i) {
    double t = i*dt;

    // true input (drifted angular velocity)
    Vector3d wt(.2, .1, 0); // this is the undrifted
    wt = wt + xts[i].bg;    // drifted angular velocity

    Vector3d at = xts[i].R.transpose()*imu.a0 + xts[i].ba;      // drifted acceleration

    Vector6d ut;
    ut << wt, at;

    // generate true
    ins.Step(xts[i+1], t, xts[i], ut, dt);

    // drifted noisy gyro reading
    Vector3d w = wt + ins.sv*Vector3d::Ones()*randn();
    // drifted noisy acc reading
    Vector3d a = at + ins.sra*Vector3d::Ones()*randn();

    Vector6d u;
    u << w, a;

    InsState x;

    timer_start(timer);
    kp.Predict(x, t, xs[i], u, dt);
    long us = timer_us(timer);
    cout << "Predict took " << us << " us." << endl;

    // this is just for printing the angles
    SO3::Instance().g2q(qt, xts[i+1].R);
    SO3::Instance().g2q(q, x.R);
    cout << "rpyt: " << qt.transpose() << endl;
    cout << "rpye: " << q.transpose() << endl;

    // true measurements of acc and mag
    at = xts[i+1].R.transpose()*imu.a0;
    //    Vector3d mt = xts[i+1].R.transpose()*sensor.m0;

    // noisy measurements of acceleration
    Vector3d za = at + xts[i+1].ba + imu.sra*Vector3d::Ones()*randn();
    //    z.normalize();

    // noisy measurements of position
    Vector3d zp = xts[i+1].p + Vector3d(gps.sxy*randn(), gps.sxy*randn(), gps.sz*randn());

    cout << "za=" << za.transpose() << endl;
    cout << "zp=" << zp.transpose() << endl;

    cout << "P-=" << x.P << endl;

    timer_start(timer);
    InsState xp;
    kci.Correct(xp, t, x, u, za);
    kcg.Correct(xs[i+1], t, xp, u, zp);
    //xs[i+1] = x;
    us = timer_us(timer);
    cout << "Update took " << us << " us." << endl;


    cout << "True attitude: " << xts[i+1].R << endl;
    cout << "Estim attitude: " << xs[i+1].R << endl;

    cout << "True biases: " << xts[i+1].bg.transpose() << endl;
    cout << "Estim biases: " << xs[i+1].bg.transpose() << endl;

    cout << "True pos: " << xts[i+1].p << endl;
    cout << "Estim pos: " << xs[i+1].p << endl;

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
