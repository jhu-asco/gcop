#include <limits>
#include "kinbody2d.h"
#include "kinbody2dmanifold.h"
#include "rn.h"
#include "se2.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Kinbody2d::Kinbody2d() : System(Kinbody2dManifold::Instance(), Rn<3>::Instance()), 
                         d(.1, .05)
{
}

double Kinbody2d::Step(Matrix3d& xb, double t, const Matrix3d& xa,
                       const Vector3d& u, double h,
                       Matrix3d *A, Matrix3d *B)
{
  SE2 &se2 = SE2::Instance();
  Matrix3d m;
  
  se2.cay(m, h*u);
  xb = xa*m;
  
  if (A) {
    se2.cay(m, -h*u);
    se2.Ad(*A, m);
  }
 
 
  if (B) {
    se2.dcay(m, -h*u);
    *B = h*m;
  }
  return 1;
}

double Kinbody2d::Step(Matrix3d &xb, double t, const Matrix3d &xa,
                       const Vector3d &u, double h,  const VectorXd &p,
                       Matrix3d *A, Matrix3d *B, Matrix<double, 3, Dynamic> *C )
{
  if (C)
    C->setZero();
  return Step(xb, t, xa, u, h, A, B);
}
