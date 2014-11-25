#include <limits>
#include "rccar1.h"
#include "rn.h"
#include "point3dmanifold.h"

using namespace gcop;
using namespace Eigen;

Rccar1::Rccar1() : Rccar(6), msteer(0.1), 
csteer(0), mdrive(0.1), cdrive(0), ktorque(1)
{
  U.bnd = true;
  U.lb[0] = -1000;
  U.ub[0] = 1000;
  U.lb[1] = -10*(M_PI/5);
  U.ub[1] = 10*(M_PI/5);
}


double Rccar1::Step(Vector4d& xb, double t, const Vector4d& xa,
                   const Vector2d& u, double h, const VectorXd *p,
                   Matrix4d *A, Matrix42d *B, Matrix4pd *C) {
  double c = cos(xa[2]);
  double s = sin(xa[2]);
  const double &v = xa[3];

  if(p)
  {
    l = (*p)[0];//Only one parameter will add wheel - body speed conversion as second parameter #TODO
    mdrive = (*p)[1];
    cdrive = (*p)[2];
    msteer = (*p)[3];
    csteer = (*p)[4];
    ktorque = (*p)[5];
    //cout<<(*p).transpose()<<endl;
  }

  double torque = ktorque*((mdrive*(bind(u[0],0)) + cdrive) - v);
  double steer = (msteer*bind(u[1],1) + csteer);
  double tansteer = tan(steer);
  //cout<<"Steer, torque, t: "<<steer<<"\t"<<torque<<"\t"<<t<<endl;
  //cout<<"torque,ktorque, mdrive, bindu, cdrive, v, t: "<<torque<<" "<<ktorque<<" "<<mdrive<<" "<<bind(u[0],0)<<" "<<cdrive<<" "<<v<<" "<<t<<endl;

  xb[0] = xa[0] + h*c*v;
  xb[1] = xa[1] + h*s*v;
  xb[2] = xa[2] + h*tansteer/l*v;
  xb[3] = v + h*torque;
  //cout<<"Xb: "<<xb.transpose()<<endl;
  //cout<<"u: "<<u.transpose()<<endl;
  
  if (A) {
    A->setIdentity();
    (*A)(0,2) = -h*s*v;
    (*A)(0,3) = h*c;
    (*A)(1,2) = h*c*v;
    (*A)(1,3) = h*s;
    (*A)(2,3) = h*tansteer/l;
  }
  
  if (B) {
    B->setZero();
    (*B)(2,1) = h*msteer/(l*v*pow(cos(steer),2));
    (*B)(3,0) = h*mdrive;
  }
  return 1;
}
