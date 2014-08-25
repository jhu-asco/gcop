#include <limits>
#include "gunicycle.h"
#include "gunicyclemanifold.h"
#include "se2.h"
#include "rn.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Gunicycle::Gunicycle() : System(GunicycleManifold::Instance()), 
                         dx(1), dy(2)
{
}


double Gunicycle::Step(M3V2d& xb, double t, const M3V2d& xa,
                       const Vector2d& u, double h, const VectorXd *p,
                       Matrix5d *A, Matrix52d *B, Matrix5Xd *C) {
  
  xb.second = xa.second + h*u;
  
  Vector3d hxib(h*xb.second[0], h*xb.second[1], 0);
  
  SE2 &se2 = SE2::Instance();
  
  Matrix3d m;
  
  se2.cay(m, hxib);
  xb.first = xa.first*m;
  
  if (A) {
    se2.cay(m, -hxib);
    Matrix3d Adm;
    se2.Ad(Adm, m);
    A->topLeftCorner<3,3>() = Adm;
    se2.dcay(m, -hxib);
    A->topRightCorner<3,2>() = h*(m.leftCols(2));
    A->bottomLeftCorner<2,3>().setZero();
    A->bottomRightCorner<2,2>().setIdentity();
  }
  
  if (B) {
    // A should've been called too
    assert(A);
    B->topRightCorner<3,2>() = h*h*(m.leftCols<2>());
    (*B)(3,0) = h; (*B)(3,1) = 0;
    (*B)(4,0) = 0; (*B)(4,1) = h;
  }
  return 1;
}
