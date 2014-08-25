#include <limits>
#include "gunicyclecost.h"
#include "gunicyclemanifold.h"
#include "se2.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

GunicycleCost::GunicycleCost(Gunicycle &sys, double tf, const M3V2d &xf) : 
  LqCost(sys, tf, xf)
{
  Q(0,0) = .01;
  Q(1,1) = .01;
  Q(2,2) = .01;
  Q(3,3) = .1;
  Q(4,4) = .1;

  Qf(0,0) = 5;
  Qf(1,1) = 5;
  Qf(2,2) = 5;
  Qf(3,3) = 1;
  Qf(4,4) = 1;

  R(0,0) = 1;
  R(1,1) = .5;
}



double GunicycleCost::L(double t, const M3V2d& x, const Vector2d& u, double h,
                        Vector5d *Lx, Matrix5d* Lxx,
                        Vector2d *Lu, Matrix2d* Luu,
                        Matrix52d *Lxu) 
{
  Vector5d xe;
  //  Vector3d ge;
  //  SE2::Instance().plog(ge, x.first);
  // xe.head(3) = ge;

  Vector3d q;
  SE2::Instance().g2q(q, x.first);
  xe.head(3) = q;
  xe.tail(2) = x.second;
  
  Matrix3d m;
  SE2::Instance().Tg(m, x.first); 

  if (t > tf - 1e-10) {
    if (Lx) {
      *Lx = Qf*xe;
      (*Lx).head(3) = m.transpose()*(*Lx).head(3);
    }
    if (Lxx) {
      *Lxx = Qf;
      (*Lxx).topLeftCorner(3,3) = m.transpose()*(*Lxx).topLeftCorner(3,3)*m;
    }
    if (Lu)
      Lu->setZero();
    if (Luu)
      Luu->setZero();
    if (Lxu)
      Lxu->setZero();

    return xe.dot(Qf*xe)/2;
    
  } else {
    if (Lx) {
      *Lx = Q*(h*xe);
      (*Lx).head(3) = h*m.transpose()*(*Lx).head(3);
    }
    if (Lxx) {
      *Lxx = h*Q;
      (*Lxx).topLeftCorner(3,3) = h*m.transpose()*Q.topLeftCorner(3,3)*m;
    }
    if (Lu)
      *Lu = h*(R*u);
    if (Luu)
      *Luu = h*R;
    if (Lxu)
      Lxu->setZero();

    return h*(xe.dot(Q*xe) + u.dot(R*u))/2;
  }
}
