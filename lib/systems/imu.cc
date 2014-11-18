#include <limits>
#include "imu.h"
#include "imumanifold.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;


Imu::Imu() : System<ImuState, 9, 3>(ImuManifold::Instance()) 
{
  sv = 3*1e-3;
  su = 3*1e-6;
  sa = 3*1e-12;
}


Imu::~Imu()
{
}

double Imu::Step(ImuState &xb, double t, const ImuState &xa, 
                 const Vector3d &u, double dt, const VectorXd *p,
                 Matrix9d *A, Matrix93d *B, Matrix9Xd *C) {
  
  SO3 &so3 = SO3::Instance();
  
  Vector3d w = u - xa.bg;             // corrected angular velocity
  
  Matrix3d dR;
  so3.exp(dR, dt*w);

  xb.R = xa.R*dR;
  xb.bg = xa.bg;
  xb.ba = xa.ba;

  Matrix3d D;
  so3.dexp(D, -dt*w);
  
  // if A is not provided, use internal A
  if (A) {
    A->setIdentity();
    A->topLeftCorner<3,3>() = dR.transpose();
    A->block<3,3>(0,3) = -dt*D;
  }
  return 0;
}

bool Imu::Noise(Matrix9d &Q, double t, const ImuState &x, const Vector3d &u, 
                double dt, const VectorXd *p) {

  double dt2 = dt*dt;
  double dt3 = dt2*dt;
  double su2 = su*su;

  Q.topLeftCorner<3,3>().diagonal().setConstant(sv*sv*dt + su2*dt3/3);
  Q.block<3,3>(0,3).diagonal().setConstant(-su2*dt2/2);
  Q.block<3,3>(3,0).diagonal().setConstant(-su2*dt2/2);
  Q.block<3,3>(3,3).diagonal().setConstant(su2*dt);
  Q.bottomRightCorner<3,3>().diagonal().setConstant(sa*sa*dt);

  return true;
}
