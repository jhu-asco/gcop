#include <limits>
#include "car.h"
#include "rn.h"

using namespace gcop;
using namespace Eigen;

Car::Car() : System(Rn<5>::Instance()), l(1)
{
}


double Car::Step(Vector5d& xb, double t, const Vector5d& xa,
                 const Vector2d& u, double h, const VectorXd *p,
                 Matrix5d *A, Matrix52d *B, Matrix5Xd *C) {
  double c = cos(xa[2]);
  double s = sin(xa[2]);
  const double &v = xa[3];
  const double &k = xa[4];
  
  xb[0] = xa[0] + h*c*v;
  xb[1] = xa[1] + h*s*v;
  xb[2] = xa[2] + h*k*v;
  xb[3] = v + h*u[0];
  xb[4] = k + h*u[1];
  
  if (A) {
    A->setIdentity();
    (*A)(0,2) = -h*s*v;
    (*A)(0,3) = h*c;
    (*A)(1,2) = h*c*v;
    (*A)(1,3) = h*s;
    (*A)(2,3) = h*k;
    (*A)(2,4) = h*v;
  }
  
  if (B) {
    B->setZero();
    (*B)(3,0) = h;
    (*B)(4,1) = h;
  }
  return 1;
}
