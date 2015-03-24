#include <limits>
#include "rccar.h"
#include "rn.h"
#include "point3dmanifold.h"

#define sgn(x) x>=0?1:-1

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

void Rccar::StateAndControlsToFlat(VectorXd &y, const Vector4d &x, const Vector2d &u)
{
	y = x.head(2);//the first two elements (x,y) are the flat outputs
}

void Rccar::FlatToStateAndControls(Vector4d &x, Vector2d &u, const std::vector<VectorXd> &y)
{
	assert(y.size() >= 3);//The flat output and derivatives upto second order should be provided
	x.head(2) = y[0];//Xand y are copied
	x(3) = y[1].norm();//yelocity is sqrt(xdot^2 + ydot^2)
	if(x(3) < 1e-6)//Very small velocity numerical error while finding u_1
	{
		u(1) = 0;
		u(2) = 0;
		x(2) = 0;
	}
	else
	{
		int signofv = sgn(x[3]);
		x(2) = atan2(signofv*y[1][1], signofv*y[1][0]);//thetadot = atan2(ydot,xdot)
		u(1) = (y[1][0]*y[2][0] + y[1][1]*y[2][1])/((this->r)*x(3));
		u(2) = (this->l/x(3))*(y[1][0]*y[2][1]-y[1][1]*y[2][0]);
	}
	//Can print x and u for debugging:
}
