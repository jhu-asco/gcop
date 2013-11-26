#ifndef GCOP_MBSLQCOST_H
#define GCOP_MBSLQCOST_H

#include "mbsmanifold.h"
#include "lqcost.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  template <int nb, int c> class MbsLqCost : public LqCost<MbsState<nb>, MBS_DIM(nb), c> {

    typedef Matrix<double, MBS_DIM(nb), 1> Vectornd;
    typedef Matrix<double, MBS_DIM(nb), MBS_DIM(nb)> Matrixnd;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, c, c> Matrixcd;
    typedef Matrix<double, MBS_DIM(nb), c> Matrixncd;
    
  public:
    
    MbsLqCost(double tf, const MbsState<nb> &xf, bool diag = true);
    
    double L(double t, const MbsState<nb> &x, const Vectorcd &u,
             Vectornd *Lx = 0, Matrixnd *Lxx = 0,
             Vectorcd *Lu = 0, Matrixcd *Luu = 0,
             Matrixncd *Lxu = 0);
  };  
  
  template <int nb, int c> MbsLqCost<c>::MbsLqCost(double tf, const MbsState<nb> &xf, bool diag) : 
    LqCost<MbsState<nb>, MBS_DIM(nb), c>(MbsManifold<nb>::Instance(), tf, xf, diag) {
   
    this->Qf.diagonal() = Matrix<double, MBS_DIM(nb), 1>::Constant(1);
    this->R.diagonal() = Matrix<double, c, 1>::Constant(.1);
  }
  
  template <int nb, int c>
    double MbsLqCost<nb, c>::L(double t, const MbsState<nb> &x, const Matrix<double, c, 1> &u,
                               Matrix<double, MBS_DIM(nb), 1> *Lx, Matrix<double, MBS_DIM(nb), MBS_DIM(nb)> *Lxx,
                               Matrix<double, c, 1> *Lu, Matrix<double, c, c> *Luu,
                               Matrix<double, MBS_DIM(nb), c> *Lxu) {
    
    Matrix<double, MBS_DIM(nb), 1> dx;
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
}


#endif
