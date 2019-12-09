#ifndef GCOP_YAWCOST_H
#define GCOP_YAWCOST_H

#include <Eigen/Dense>
#include <iostream>
#include "cost.h"
#include "so3.h"

namespace gcop {

  using namespace Eigen;
  using namespace std;
  
  template <typename T, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic> class YawCost : public Cost<T, _nx, _nu, _np> {
  public:
  
  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;

  typedef Matrix<double, _nx, _nx> Matrixnd;
  typedef Matrix<double, _nx, _nu> Matrixncd;
  typedef Matrix<double, _nu, _nx> Matrixcnd;
  typedef Matrix<double, _nu, _nu> Matrixcd;
  
  //  typedef Matrix<double, Dynamic, 1> Vectormd;
  typedef Matrix<double, _np, _np> Matrixmd;
  typedef Matrix<double, _nx, _np> Matrixnmd;
  typedef Matrix<double, _np, _nx> Matrixmnd;
  
  /**
   * create a cost interface
   * @param X the state manifold which is used to perform addition/subtraction of states
   * @param tf time horizon: when the cost function L is called it will internally check whether its argument t is equation to tf and return the terminal cost
   */
  YawCost(System<T, _nx, _nu, _np> &sys, double tf, Vector3d targ);
  
  /**
   * Cost function L
   * @param t time
   * @param x state
   * @param u control
   * @param h time-step
   * @param Lx derivative wrt x
   * @param Lxx derivative wrt x,x
   * @param Lu derivative wrt u
   * @param Luu derivative wrt u,u
   * @param Lxu derivative wrt x,u     
   */
  virtual double L(double t, const T &x, const Vectorcd &u, double h,
                   const Vectormd *p = 0,
                   Vectornd *Lx = 0, Matrixnd* Lxx = 0,
                   Vectorcd *Lu = 0, Matrixcd* Luu = 0,
                   Matrixncd *Lxu = 0,
                   Vectormd *Lp = 0, Matrixmd *Lpp = 0,
                   Matrixmnd *Lpx = 0);  

  Vector3d target;
  double gain;
  
  };


  template <typename T, int _nx, int _nu, int _np> 
    YawCost<T, _nx, _nu, _np>::YawCost(System<T, _nx, _nu, _np> &sys, 
                                           double tf, Vector3d t) : Cost<T, _nx, _nu, _np>(sys, tf) {
    target = t;
    gain = 1;
  }  
  
  template <typename T, int _nx, int _nu, int _np> 
    double YawCost<T, _nx, _nu, _np>::L(double t, const T& x, const Vectorcd& u, double h,
                                          const Vectormd *p,
                                          Matrix<double, _nx, 1> *Lx, Matrix<double, _nx, _nx>* Lxx,
                                          Matrix<double, _nu, 1> *Lu, Matrix<double, _nu, _nu>* Luu,
                                          Matrix<double, _nx, _nu> *Lxu,
                                          Vectormd *Lp, Matrixmd *Lpp,
                                          Matrixmnd *Lpx) {
    cout << "Yaw L Called" << endl;
    Vector3d e(0,0,1);
    Vector3d df = target - x.p;
    cout << "target: " << target << endl;
    cout << "pos: " << x.p << endl;
    cout << "df: " << df << endl;
    df(2) = 0;
    cout << "df: " << df << endl;
    double dnorm = df.norm();
    cout << "df Norm: " << dnorm << endl;
    double dx = df(0);
    double dy = df(1);
    df = df/df.norm();
    double c = gain*(1 - (x.R*e).dot(df));
    if (dnorm < 1e-10)
    { return c; }
    Matrix3d ehat;
    SO3::Instance().hat(ehat, e);
    Matrix3d ddfdp;
    ddfdp << ((-1/dnorm) + dx*dx/(dnorm*dnorm*dnorm)), (dx*dy)/(dnorm*dnorm*dnorm), 0,
             (dx*dy)/(dnorm*dnorm*dnorm), ((-1/dnorm) + dy*dy/(dnorm*dnorm*dnorm)), 0,
             0, 0, 0;
    MatrixXd d2dfdp2(9,3);
    d2dfdp2 << (Matrix2d()<<(3*dx/(dnorm*dnorm*dnorm) - 3*(dx*dx*dx)/(dnorm*dnorm*dnorm*dnorm*dnorm)),
                            (dy/(dnorm*dnorm*dnorm) - 3*(dx*dx*dy)/(dnorm*dnorm*dnorm*dnorm*dnorm)),
                            (dy/(dnorm*dnorm*dnorm) - 3*(dx*dx*dy)/(dnorm*dnorm*dnorm*dnorm*dnorm)),
                            (dx/(dnorm*dnorm*dnorm) - 3*(dx*dy*dy)/(dnorm*dnorm*dnorm*dnorm*dnorm))).finished(),
               MatrixXd::Zero(2,1), MatrixXd::Zero(1,3),
               (Matrix2d()<<(dy/(dnorm*dnorm*dnorm) - 3*(dx*dx*dy)/(dnorm*dnorm*dnorm*dnorm*dnorm)),
                             (dx/(dnorm*dnorm*dnorm) - 3*(dx*dy*dy)/(dnorm*dnorm*dnorm*dnorm*dnorm)),
                             (dx/(dnorm*dnorm*dnorm) - 3*(dx*dx*dx)/(dnorm*dnorm*dnorm*dnorm*dnorm)),
                             (3*dy/(dnorm*dnorm*dnorm) - 3*(dx*dx*dy)/(dnorm*dnorm*dnorm*dnorm*dnorm))).finished(),
               MatrixXd::Zero(2,1), MatrixXd::Zero(4,3);
    MatrixXd Rstack(3,9);
    Rstack << (x.R*e).transpose(), MatrixXd::Zero(1,6),
              MatrixXd::Zero(1,3), (x.R*e).transpose(), MatrixXd::Zero(1,3),
              MatrixXd::Zero(1,9);
    if (Lx){
    cout << "Yaw Lx Called" << endl;
      Lx->setZero();
      Vector3d dLdR = c*(ehat.transpose()*x.R.transpose()*df);
      Vector3d dLdp = -c*(ddfdp.transpose()*x.R*e);
      Lx->head(3) = dLdR;
      Lx->segment(3,3) = dLdp;
    }
    if (Lxx){
    cout << "Yaw Lxx Called" << endl;
      Lxx->setZero();
      Matrix3d RTdfhat;
      SO3::Instance().hat(RTdfhat,x.R.transpose()*df);
      Matrix3d d2LdR2 = -c*ehat*RTdfhat;
      Matrix3d d2LdRdp = c*(ehat.transpose()*x.R.transpose()*ddfdp);
      Matrix3d d2LdpdR = c*(ddfdp.transpose()*x.R*ehat);
      Matrix3d d2Ldp2 = -c*(Rstack*d2dfdp2);
      Lxx->topLeftCorner(3,3) = d2LdR2;
      Lxx->block(0,3,3,3) = d2LdRdp;
      Lxx->block(3,0,3,3) = d2LdpdR;
      Lxx->block(3,3,3,3) = d2Ldp2;
    }
    if (Lu)
      Lu->setZero();
    if (Luu)
      Luu->setZero();
    if (Lxu)
      Lxu->setZero();
    if (Lp)
      Lp->setZero();
    if (Lpp)
      Lpp->setZero();
    if (Lpx)
      Lpx->setZero();

    cout << "Yaw L Returned" << endl;
    return c;
  }
}

#endif
