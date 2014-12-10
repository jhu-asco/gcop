#include <limits>
#include "rccar.h"
#include "rn.h"
#include "point3dmanifold.h"

using namespace gcop;
using namespace Eigen;

Rccar::Rccar(int np) : System(Rn<4>::Instance(),2,np), l(.3), r(0.5), h(0.01)
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

  if(p)
  {
    l = (*p)[0];//length of car
    r = (*p)[1];//radius of the car
  }
 const double &force = r*u[0];
 const double &tansteer = u[1];
 int nofsteps = floor(h/(this->h));
 Vector4d xtemp = xa;
 if (nofsteps > 0)
 {
   const double &v = xtemp[3];
   for(int count = 0;count < nofsteps; count++)
   {
     xtemp[0] = xtemp[0] + (this->h)*cos(xtemp[2])*v;
     xtemp[1] = xtemp[1] + (this->h)*sin(xtemp[2])*v;
     xtemp[2] = xtemp[2] + (this->h)*tansteer/l*v;
     xtemp[3] = v + (this->h)*force;
   }
 }
 double remh = h - nofsteps*(this->h);//Left remainder

 xb[0] = xtemp[0] + remh*cos(xtemp[2])*xtemp[3];
 xb[1] = xtemp[1] + remh*sin(xtemp[2])*xtemp[3];
 xb[2] = xtemp[2] + remh*tansteer/l*xtemp[3];
 xb[3] = xtemp[3] + remh*force;

  const double &v = xa[3];
  double c = cos(xa[2]);
  double s = sin(xa[2]);
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
