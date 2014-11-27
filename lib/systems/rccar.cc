#include <limits>
#include "rccar.h"
#include "rn.h"
#include "point3dmanifold.h"

using namespace gcop;
using namespace Eigen;

Rccar::Rccar(int np) : System(Rn<4>::Instance(),2,np), l(.3), r(0.5)
{
  U.bnd = true;
  U.lb[0] = -100;
  U.ub[0] = 100;
  U.lb[1] = tan(-M_PI/5);
  U.ub[1] = tan(M_PI/5);
}


double Rccar::Step(Vector4d& xb, double t, const Vector4d& xa,
                   const Vector2d& u, double h, const VectorXd *p,
                   Matrix4d *A, Matrix42d *B, Matrix4pd *C) {
  double c = cos(xa[2]);
  double s = sin(xa[2]);
  const double &v = xa[3];

  if(p)
  {
    l = (*p)[0];//length of car
    r = (*p)[1];//radius of the car
  }
 const double &force = r*u[0];
 const double &tansteer = u[1];
  xb[0] = xa[0] + h*c*v;
  xb[1] = xa[1] + h*s*v;
  xb[2] = xa[2] + h*tansteer/l*v;
  xb[3] = v + h*force;
  
  if (A) {
    A->setIdentity();
    (*A)(0,2) = -h*s*v;
    (*A)(0,3) = h*c;
    (*A)(1,2) = h*c*v;
    (*A)(1,3) = h*s;
    (*A)(2,3) = h*u[1]/l;
  }
  
  if (B) {
    B->setZero();
    (*B)(2,1) = h/l*v;
    (*B)(3,0) = h*r;
  }
  return 1;
}
