#include <limits>
#include "gcar.h"
#include "gcarmanifold.h"
#include "se2.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Gcar::Gcar() : System(GcarManifold::Instance()), 
               l(.5), r(1)
{
  U.bnd = true;
  U.lb[0] = -100;
  U.ub[0] = 100;
  U.lb[1] = tan(-M_PI/5);
  U.ub[1] = tan(M_PI/5);
}


double Gcar::Step(M3V1d& xb, double t, const M3V1d& xa,
                  const Vector2d& u, double h, const VectorXd *p,
                  Matrix4d *A, Matrix42d *B, Matrix4Xd *C) {
  
  // update forward velocity
  xb.second = xa.second + h*r*u[0];
  
  // the body-fixed velocity vector multiplied by delta-t
  Vector3d hxib(h*u[1]/l*xb.second, h*xb.second, 0);
  
  SE2 &se2 = SE2::Instance();
    
  // update configuration on SE(2)
  Matrix3d m;
  se2.cay(m, hxib);
  xb.first = xa.first*m;

  /*
  if (A) {
    
    se2.cay(m, -hxib);
    Matrix3d Adm;
    se2.Ad(Adm, m);
    A->topLeftCorner<3,3>() = Adm;
    se2.dcay(m, -hxib);
    A->topRightCorner<3,1>() = m*Vector3d(h*u[1]/l, h, 0);
    A->bottomLeftCorner<1,3>().setZero();
    (*A)(3,3) = 1;
  }
  
  if (B) {
    // A should've been called too
    assert(A);
    Matrix<double, 3, 2> F;
    F(0,0) = h*h*u[1]/l*r;
    F(0,1) = h/l*xb.second;
    F(1,0) = h*h*r;
    F(1,1) = F(2,0) = F(2,1) = 0;
    B->topRightCorner<3,2>() = m*F;
    (*B)(3,0) = r*h; (*B)(3,1) = 0;
  }
  */
  return 1;
}
