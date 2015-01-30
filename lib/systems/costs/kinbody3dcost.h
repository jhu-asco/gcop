#ifndef GCOP_KINBODY3DCOST_H
#define GCOP_KINBODY3DCOST_H

#include "kinbody3d.h"
#include "kinbody3dtrack.h"
#include "lqcost.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  class Kinbody3dCost : public LqCost<Matrix4d, 6, 6> {
    
    typedef Matrix<double, 6, 1> Vector6d;
    typedef Matrix<double, 6, 6> Matrix6d;
    
  public:
    
    Kinbody3dCost(Kinbody3d &sys, double tf, const Matrix4d &xf, bool diag = true);
    
    /*    double L(double t, const Body3dState &x, const Vectorcd &u,
             Vector12d *Lx = 0, Matrix12d *Lxx = 0,
             Vectorcd *Lu = 0, Matrixcd *Luu = 0,
             Matrix12xcd *Lxu = 0);
    */
    Kinbody3dTrack *track;
    double ko;
  };  
  
  Kinbody3dCost::Kinbody3dCost(Kinbody3d &sys, double tf, const Matrix4d &xf, bool diag) : 
    LqCost<Matrix4d, 6, 6>(sys, tf, xf, diag) {

    /*
      Q(0,0) = .05;
      Q(1,1) = .05;
      Q(2,2) = .05;
      Q(3,3) = .1;
      Q(4,4) = .1;
      Q(5,5) = .1;
      Q(6,6) = 0;
      Q(7,7) = 0;
      Q(8,8) = 0;
      Q(9,9) = 0;
      Q(10,10) = 0;
      Q(11,11) = 0;
    */
    
    this->Qf(0,0) = .5;
    this->Qf(1,1) = .5;
    this->Qf(2,2) = .5;
    this->Qf(3,3) = 5;
    this->Qf(4,4) = 5;
    this->Qf(5,5) = 5;
    
    this->R.diagonal() = Matrix<double, 6, 1>::Constant(.1);
  }
  /*  
template <int c>
  double Body3dCost<c>::L(double t, const Body3dState &x, const Matrix<double, c, 1> &u,
                          Vector12d *Lx, Matrix12d *Lxx,
                          Matrix<double, c, 1> *Lu, Matrix<double, c, c> *Luu,
                          Matrix<double, 12, c> *Lxu)
{
  Vector12d dx;
  this->X.Lift(dx, this->xf, x);

  
  // check if final state
  if (t > this->tf - 1e-10) {
    if (Lx) {
      if (this->diag)
        *Lx = this->Qf.diagonal().cwiseProduct(dx);
      else
        *Lx = this->Qf*dx;
      
      // add dcayinv if this->Q(1:3,1:3) != a*Id
      //      (*Lx).head(3) = Rt*(*Lx).head<3>();
    }
    if (Lxx) {
      *Lxx = this->Qf;
      //      (*Lxx).topLeftCorner(3,3) = Rt*(*Lxx).topLeftCorner<3,3>()*R;
    }

    if (Lu)
      Lu->setZero();
    if (Luu)
      Luu->setZero();
    if (Lxu)
      Lxu->setZero();

    if (this->diag)
      return dx.dot(this->Qf.diagonal().cwiseProduct(dx))/2;
    else
      return dx.dot(this->Qf*dx)/2;
    
  } else {
    if (Lx) {
      if (this->diag)
        *Lx = this->Q.diagonal().cwiseProduct(dx);
      else
        *Lx = this->Q*dx;
      //      (*Lx).head<3>() = Rat*(*Lx).head<3>();
    }
    if (Lxx) {
      *Lxx = this->Q;
      //      (*Lxx).topLeftCorner<3,3>() = Rt*this->Q.topLeftCorner<3,3>()*R;
    }
    if (Lu)
      if (this->diag)
        *Lu = this->R.diagonal().cwiseProduct(u);
      else
        *Lu = this->R*u;

    if (Luu)
      *Luu = this->R;
      
    if (Lxu)
      Lxu->setZero();

    if (this->diag)
      return (dx.dot(this->Q.diagonal().cwiseProduct(dx)) + u.dot(this->R.diagonal().cwiseProduct(u)))/2;
    else
      return (dx.dot(this->Q*dx) + u.dot(this->R*u))/2;
  }
}
  */
}


#endif
