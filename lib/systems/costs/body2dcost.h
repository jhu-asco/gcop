#ifndef GCOP_BODY2DCOST_H
#define GCOP_BODY2DCOST_H

#include "lqcost.h"
#include <limits>
#include "body2dtrack.h"
#include "body2d.h"
#include "se2.h"
#include <iostream>
#include <assert.h>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 6, 3> Matrix63d;
  typedef pair<Matrix3d, Vector3d> Body2dState;
  
  template<int c = 3>
    class Body2dCost : public LqCost<Body2dState, 6, c, Dynamic, 9> {
  public:
  typedef Matrix<double, c, 1> Vectorcd;
  typedef Matrix<double, 6, c> Matrix6cd;
  typedef Matrix<double, 3, c> Matrix3cd;
  typedef Matrix<double, c, c> Matrixcd;
  
  Body2dCost(Body2d<c> &sys, double tf, const Body2dState &xf);
    
  double L(double t, const Body2dState &x, const Vectorcd &u, double h,
           Vector6d *Lx = 0, Matrix6d *Lxx = 0,
           Vectorcd *Lu = 0, Matrixcd *Luu = 0,
           Matrix6cd *Lxu = 0);
  
  
  Body2dTrack *track;
  double ko;  ///< obstacle avoidance gain
  
  };
  
 
 template<int c>
   Body2dCost<c>::Body2dCost(Body2d<c> &sys, double tf, const Body2dState &xf) : 
   LqCost<Body2dState, 6, c, Dynamic, 9>(sys, tf, xf), track(0), ko(.3)
   {
     this->Q(0,0) = .01;
     this->Q(1,1) = .01;
     this->Q(2,2) = .01;
     this->Q(3,3) = .1;
     this->Q(4,4) = .1;
     this->Q(5,5) = .1;
     
     this->Qf(0,0) = 5;
     this->Qf(1,1) = 5;
     this->Qf(2,2) = 5;
     this->Qf(3,3) = 1;
     this->Qf(4,4) = 1;
     this->Qf(5,5) = 1;

     this->R.setIdentity();
   }
 


 template<int c>
   double Body2dCost<c>::L(double t, const Body2dState& x, const Vectorcd& u, double h,
                           Vector6d *Lx, Matrix6d* Lxx,
                           Vectorcd *Lu, Matrixcd* Luu,
                           Matrix6cd *Lxu) 
   {
     Vector6d xe;
     //  Vector3d ge;
     //  SE2::Instance().plog(ge, x.fthis->irst);
     // xe.head(3) = ge;
     
     Vector3d q;
     Matrix3d gi;
     
     
     Matrix2d P;
     bool inside = false;
     if (track && 0) {    
       Vector2d p = x.first.topRightCorner<2,1>(); //current position
       
       //    Vector2d s = track->w/2*this->xf->first.block<2,1>(0,1); // from center to a
       Vector2d s = track->w/2*this->xf->first.block(0,1,2,1); // from center to a
       Vector2d a(this->xf->first(0,2) + s(0), this->xf->first(1,2) + s(1));
       Vector2d b(this->xf->first(0,2) - s(0), this->xf->first(1,2) - s(1));
       Vector2d v = b - a;
       double vn = v.norm();
       v = v/vn;
       
       Matrix3d gp = this->xf->first; // projected position
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
       SE2::Instance().inv(gi, this->xf->first);
     }
     
     Matrix3d dg = gi*x.first;
     SE2::Instance().g2q(q, dg);
     xe.head(3) = q;
     xe.tail(3) = x.second;
     
     Matrix3d m;
     SE2::Instance().Tg(m, dg);
     
     if (t > this->tf - 1e-10) {
       if (Lx) {
         *Lx = this->Qf*xe;
         if (track && inside)
           Lx->segment<2>(1) = P*Lx->segment<2>(1);
         (*Lx).head(3) = m.transpose()*(*Lx).head(3);
       }
       if (Lxx) {
         *Lxx = this->Qf;
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
       
       return xe.dot(this->Qf*xe)/2;
       
     } else {
       if (Lx) {
         *Lx = this->Q*(h*xe);
         (*Lx).head(3) = h*m.transpose()*(*Lx).head(3);
       }
       if (Lxx) {
         *Lxx = h*this->Q;
         (*Lxx).topLeftCorner(3,3) = h*m.transpose()*this->Q.topLeftCorner(3,3)*m;
       }
       if (Lu)
         *Lu = h*(this->R*u);
       if (Luu)
         *Luu = h*this->R;
       if (Lxu)
         Lxu->setZero();
       
       double e = h*(xe.dot(this->Q*xe) + u.dot(this->R*u))/2;
       
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
             
             double w = 1/d - 1/dmax;
             eo += ko/2*w*w;
             double d2 = d*d;
             double d3 = d2*d;
             double d4 = d3*d;
             const Matrix2d &R = x.first.topLeftCorner<2,2>();
             if (Lx)
               Lx->tail<2>() += R.transpose()*(-ko*w*v/d3);
             if (Lxx)
               Lxx->bottomRightCorner<2,2>() += R.transpose()*( (ko*(1/d+3/d2-1/dmax-2/(dmax*d))/d4*v)*v.transpose() - (ko*w/d3)*Matrix2d::Identity())*R;
           }
         } else {
           for (int j = 0; j < track->ls.size(); ++j) {
             Vector2d v = p - track->ls[closest ? l : j];
             double d = v.norm();
             if (d > dmax)
               continue;
             
             double w = 1/d - 1/dmax;
             eo += ko/2*w*w;
             double d2 = d*d;
             double d3 = d2*d;
             double d4 = d3*d;
             const Matrix2d &R = x.first.topLeftCorner<2,2>();
             if (Lx)
               Lx->tail<2>() += R.transpose()*(-ko*w*v/d3);
             if (Lxx)
               Lxx->bottomRightCorner<2,2>() += R.transpose()*( (ko*(1/d+3/d2-1/dmax-2/(dmax*d))/d4*v)*v.transpose() - (ko*w/d3)*Matrix2d::Identity())*R;
             
             if (closest)
               break;
           }
         }
         cout << "obstacle cost=" << eo << " t=" << t << " tf=" << this->tf << endl;
         e += eo;      
       }
       return e;
     }   
   } 
}

#endif
