#include <limits>
#include "kinbody3d.h"
#include "kinbody3dmanifold.h"
#include "rn.h"
#include "se3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Kinbody3d::Kinbody3d() : System(Kinbody3dManifold::Instance()), 
  d(.1, .05, .05)
{
}

double Kinbody3d::Step(Matrix4d &xb, double t, const Matrix4d &xa,
                       const Vector6d &u, double h,  const VectorXd *p,
                       Matrix6d *A, Matrix6d *B, Matrix<double, 6, Dynamic> *C )
{
  SE3 &se3 = SE3::Instance();
  Matrix4d m;
  
  se3.cay(m, h*u);
  xb = xa*m;
  
  if (A) {
    se3.cay(m, -h*u);
    se3.Ad(*A, m);
  }
  
  if (B) {
    Matrix6d mb;
    se3.dcay(mb, -h*u);
    *B = h*mb;
  }

  if (C)
    C->setZero();
}
