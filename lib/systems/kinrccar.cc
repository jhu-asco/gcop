#include <limits>
#include "kinrccar.h"
#include "se3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

KinRccar::KinRccar() : Kinbody3d<2>()
{
  d << 0.3, 0.25, 0.2;
}

double KinRccar::Step(Matrix4d &xb, double t, const Matrix4d &xa,
                       const Vector2d &u, double h,  const VectorXd *p,
                       Matrix6d *A, Matrix62d *B, Matrix<double, 6, Dynamic> *C )
{
  SE3 &se3 = SE3::Instance();
  Matrix4d m;
  
  Vector6d u_kin;
  u_kin << 0, 0, u(0)*tan(u(1))/d(0), u(0), 0, 0; 

  se3.cay(m, h*u_kin);
  xb = xa*m;
  
  if (A) {
    A->setZero();
    se3.cay(m, -h*u_kin);
    se3.Ad(*A, m);
  }
  
  if (B) {
    B->setZero();
    Matrix6d mb;
    Matrix62d jac;
    jac << 0, 0,
           0, 0,
           tan(u(1))/d(0), 1./(d(0)*cos(u(1))*cos(u(1))), 
           1, 0,
           0, 0,
           0, 0;
    se3.dcay(mb, -h*u_kin);
    *B = h*mb*jac;
  }

  if (C)
  {
    C->setZero();
  }
}
