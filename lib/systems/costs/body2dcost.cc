#include <limits>
#include "body2dcost.h"
#include "body2dmanifold.h"
#include "se2.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Body2dCost::Body2dCost(Body2d &sys, double tf, const Body2dState &xf) : 
  LqCost(sys, tf, xf), track(0), ko(.3)
{
  Q(0,0) = .01;
  Q(1,1) = .01;
  Q(2,2) = .01;
  Q(3,3) = .1;
  Q(4,4) = .1;
  Q(5,5) = .1;

  Qf(0,0) = 5;
  Qf(1,1) = 5;
  Qf(2,2) = 5;
  Qf(3,3) = 1;
  Qf(4,4) = 1;
  Qf(5,5) = 1;

  R(0,0) = 1;
  R(1,1) = .5;
  R(2,2) = .5;
}


double Body2dCost::L(double t, const Body2dState& x, const Vector3d& u, double h,
                     Vector6d *Lx, Matrix6d* Lxx,
                     Vector3d *Lu, Matrix3d* Luu,
                     Matrix63d *Lxu) 
{
  Vector6d xe;
  //  Vector3d ge;
  //  SE2::Instance().plog(ge, x.first);
  // xe.head(3) = ge;

  Vector3d q;
  Matrix3d gi;


  Matrix2d P;
  bool inside = false;
  if (track && 0) {    
    Vector2d p = x.first.topRightCorner<2,1>(); //current position
    
    Vector2d s = track->w/2*xf->first.block<2,1>(0,1); // from center to a
    Vector2d a(xf->first(0,2) + s(0), xf->first(1,2) + s(1));
    Vector2d b(xf->first(0,2) - s(0), xf->first(1,2) - s(1));
    Vector2d v = b - a;
    double vn = v.norm();
    v = v/vn;
    
    Matrix3d gp = xf->first; // projected position
    double r = v.transpose()*(p - a);
    if (r < 0) {
      gp.topRightCorner<2,1>() = a;
    } else
      if (r > vn) {
          gp.topRightCorner<2,1>() = b;
        } else {
          inside = true;
          gp.topRightCorner<2,1>() = a + v*r;
          SE2::Instance().inv(gi, gp);
          P = Matrix2d::Identity() - v*v.transpose();
        }
  } else {
    SE2::Instance().inv(gi, xf->first);
  }
  
  Matrix3d dg = gi*x.first;
  SE2::Instance().g2q(q, dg);
  xe.head(3) = q;
  xe.tail(3) = x.second;
  
  Matrix3d m;
  SE2::Instance().Tg(m, dg);

  if (t > tf - 1e-10) {
    if (Lx) {
        *Lx = Qf*xe;
        if (track && inside)
          Lx->segment<2>(1) = P*Lx->segment<2>(1);
        (*Lx).head(3) = m.transpose()*(*Lx).head(3);
    }
    if (Lxx) {
      *Lxx = Qf;
      if (track && inside)
        Lxx->block<2,2>(1,1) = P*Lxx->block<2,2>(1,1)*P;      
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

    double e = h*(xe.dot(Q*xe) + u.dot(R*u))/2;
  
    if (track) {
      double dmax = 1;
     
      double eo = 0;
      Vector2d p(x.first(0,2), x.first(1,2));

      bool closest = true;     
      int l = 0;

      double dmin = 1e16;
      if (closest) {
        for (int j = 0; j < track->ls.size(); ++j) {
          Vector2d v = p - track->ls[j];
          double d = v.norm();
          if (d < dmin) {
            dmin = d;
            l = j;
          }
        }
        
        if (dmin < dmax)  {
          Vector2d v = p - track->ls[l];
          double d = v.norm();
          
          double c = 1/d - 1/dmax;
          eo += ko/2*c*c;
          double d2 = d*d;
          double d3 = d2*d;
          double d4 = d3*d;
          const Matrix2d &R = x.first.topLeftCorner<2,2>();
          if (Lx)
            Lx->tail<2>() += R.transpose()*(-ko*c*v/d3);
          if (Lxx)
            Lxx->bottomRightCorner<2,2>() += R.transpose()*( (ko*(1/d+3/d2-1/dmax-2/(dmax*d))/d4*v)*v.transpose() - (ko*c/d3)*Matrix2d::Identity())*R;
        }
      } else {
        for (int j = 0; j < track->ls.size(); ++j) {
          Vector2d v = p - track->ls[closest ? l : j];
          double d = v.norm();
          if (d > dmax)
            continue;
          
          double c = 1/d - 1/dmax;
          eo += ko/2*c*c;
          double d2 = d*d;
          double d3 = d2*d;
          double d4 = d3*d;
          const Matrix2d &R = x.first.topLeftCorner<2,2>();
          if (Lx)
            Lx->tail<2>() += R.transpose()*(-ko*c*v/d3);
          if (Lxx)
            Lxx->bottomRightCorner<2,2>() += R.transpose()*( (ko*(1/d+3/d2-1/dmax-2/(dmax*d))/d4*v)*v.transpose() - (ko*c/d3)*Matrix2d::Identity())*R;
          
          if (closest)
            break;
        }
      }
      cout << "obstacle cost=" << eo << " t=" << t << " tf=" << tf << endl;
      e += eo;      
    }
    return e;
  }   
}
