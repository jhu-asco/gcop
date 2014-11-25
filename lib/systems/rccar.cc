#include <limits>
#include "rccar.h"
#include "rn.h"
#include "point3dmanifold.h"

using namespace gcop;
using namespace Eigen;

Rccar::Rccar(int np) : 
#ifdef RCCAR_RADIO_INPUTS
System(Rn<4>::Instance(),2,np), l(.3),msteer(0.1), 
csteer(0), mdrive(0.1), cdrive(0), ktorque(1)

#else
System(Rn<4>::Instance(),2,np), l(.3)
#endif
{
  U.bnd = true;
#ifdef  RCCAR_RADIO_INPUTS
  U.lb[0] = -1000;
  U.ub[0] = 1000;
  U.lb[1] = -10*(M_PI/5);
  U.ub[1] = 10*(M_PI/5);
#else
  U.lb[0] = -100;
  U.ub[0] = 100;
  U.lb[1] = tan(-M_PI/5);
  U.ub[1] = tan(M_PI/5);
#endif
}


double Rccar::Step(Vector4d& xb, double t, const Vector4d& xa,
                   const Vector2d& u, double h, const VectorXd *p,
                   Matrix4d *A, Matrix42d *B, Matrix4pd *C) {
  double c = cos(xa[2]);
  double s = sin(xa[2]);
  const double &v = xa[3];

  if(p)
  {
    l = (*p)[0];//Only one parameter will add wheel - body speed conversion as second parameter #TODO
#ifdef RCCAR_RADIO_INPUTS
    mdrive = (*p)[1];
    msteer = (*p)[2];
    cdrive = (*p)[3];
    csteer = (*p)[4];
    ktorque = (*p)[5];
#endif
  }
#ifdef RCCAR_RADIO_INPUTS
 double torque = ktorque*((mdrive*(bind(u[0],0)) + cdrive) - v);
 double steer = (msteer*bind(u[1],1) + csteer);
 double tansteer = tan(steer);
#else
 const double &torque = u[0];
 const double &tansteer = u[1];
#endif
  xb[0] = xa[0] + h*c*v;
  xb[1] = xa[1] + h*s*v;
  xb[2] = xa[2] + h*tansteer/l*v;
  xb[3] = v + h*torque;
  
  if (A) {
    A->setIdentity();
    (*A)(0,2) = -h*s*v;
    (*A)(0,3) = h*c;
    (*A)(1,2) = h*c*v;
    (*A)(1,3) = h*s;
#ifdef RCCAR_RADIO_INPUTS
    (*A)(2,3) = h*tansteer/l;
#else
    (*A)(2,3) = h*u[1]/l;
#endif
  }
  
  if (B) {
    B->setZero();
#ifdef RCCAR_RADIO_INPUTS
    (*B)(2,1) = h*msteer/(l*v*pow(cos(steer),2));
    (*B)(3,0) = h*mdrive;
#else
    (*B)(2,1) = h/l*v;
    (*B)(3,0) = h;
#endif
  }
  return 1;
}
