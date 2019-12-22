#ifndef GCOP_YAWCOST_H
#define GCOP_YAWCOST_H

#include <Eigen/Dense>
#include <iostream>
#include "cost.h"
#include "so3.h"
#include <iomanip>
#include <iostream>

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
  YawCost(System<T, _nx, _nu, _np> &sys, double tf, const Vector3d &targ);
  
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
  const Vector3d e = Vector3d::UnitX();
  double gain;
  //bool useFinite;

  double simple_c(const T &x);
  void pop_Lx(const T &x, double h, Matrix<double, _nx, 1> *Lx);
  };

  template <typename T, int _nx, int _nu, int _np> 
    double YawCost<T, _nx, _nu, _np>::simple_c(const T& x) {

    Vector3d df = target - x.p;
    //df(2) = 0;
    double dnorm = df.norm();
    if (dnorm < 1e-10) {
      return 0;
    }
    df = df/df.norm();
    return (gain*(1 - (x.R*e).dot(df)));
  }
  template <typename T, int _nx, int _nu, int _np>
    void YawCost<T, _nx, _nu, _np>::pop_Lx(const T& x, double h, Matrix<double, _nx, 1> *Lx){
      Vector3d df = target - x.p;
      //df(2) = 0;
      double dnorm = df.norm();
      if (dnorm < 1e-10) {
        return;
      }
      double dx = df(0);
      double dy = df(1);
      double dz = df(2);
      df = df/df.norm();
      Matrix3d ehat;
      SO3::Instance().hat(ehat, e);
      Matrix3d ddfdp;
      ddfdp << (dy*dy + dz*dz), (-dx*dy), (-dx*dz),
               (-dx*dy), (dx*dx + dz*dz), (-dy*dz),
               (-dx*dz), (-dy*dz), (dx*dx + dy*dy);
      ddfdp = ddfdp/(dnorm*dnorm*dnorm);
      Vector3d dLdR = gain*(ehat.transpose()*x.R.transpose()*df);
      Vector3d dLdp = gain*(ddfdp.transpose()*x.R*e);
      Lx->setZero();
      Lx->head(3) = dLdR;
      Lx->segment(3,3) = dLdp;
      /*double c = simple_c(x);
      double eps = 1e-10;
      Vectornd dx;
      for (int i = 0; i < Cost<T, _nx, _nu, _np>::sys.X.n; ++i ) {
        //Set dx to eps in the the i-th dim
        dx.setZero();
        dx[i] = eps;
        //Calculate using dx
        T xtemp;
        Cost<T, _nx, _nu, _np>::sys.X.Retract(xtemp,x,dx);//xtemp = xa+dx
        Cost<T, _nx, _nu, _np>::sys.Rec(xtemp,h);//Reconstruct
        double dpos = simple_c(xtemp) - c;//dpos = c(xtemp) - c(x)
        //Calculate using -dx
        Cost<T, _nx, _nu, _np>::sys.X.Retract(xtemp,x,-dx);//xtemp = xa+dx
        Cost<T, _nx, _nu, _np>::sys.Rec(xtemp,h);//Reconstruct
        double dneg = simple_c(xtemp) - c;//dneg = c(xtemp) - c(x)
        //Set Lx[i]
        (*Lx)[i] = (dpos - dneg)/(2*eps);
      }*/
  }

  template <typename T, int _nx, int _nu, int _np> 
    YawCost<T, _nx, _nu, _np>::YawCost(System<T, _nx, _nu, _np> &sys, 
                                       double tf, const Vector3d &targ) : 
       Cost<T, _nx, _nu, _np>(sys, tf), target(targ), gain(1) {
  }  
  
  template <typename T, int _nx, int _nu, int _np> 
    double YawCost<T, _nx, _nu, _np>::L(double t, const T& x, const Vectorcd& u, double h,
                                          const Vectormd *p,
                                          Matrix<double, _nx, 1> *Lx, Matrix<double, _nx, _nx>* Lxx,
                                          Matrix<double, _nu, 1> *Lu, Matrix<double, _nu, _nu>* Luu,
                                          Matrix<double, _nx, _nu> *Lxu,
                                          Vectormd *Lp, Matrixmd *Lpp,
                                          Matrixmnd *Lpx) {
    Vector3d df = target - x.p;
    //df(2) = 0;
    double dnorm = df.norm();
    if (dnorm < 1e-10) {
      return 0;
    }
    double dx = df(0);
    double dy = df(1);
    double dz = df(2);
    df = df/df.norm();
    double c = simple_c(x);
    Matrix3d ehat;
    SO3::Instance().hat(ehat, e);
    Matrix3d ddfdp;
    ddfdp << (dy*dy + dz*dz), (-dx*dy), (-dx*dz),
             (-dx*dy), (dx*dx + dz*dz), (-dy*dz),
             (-dx*dz), (-dy*dz), (dx*dx + dy*dy);
    ddfdp = ddfdp/(dnorm*dnorm*dnorm);
    MatrixXd d2dfdp2(9,3);
    d2dfdp2 << -(3*dx*(dy*dy + dz*dz)), -(dy*(-2*dx*dx + dy*dy + dz*dz)), -(dz*(-2*dx*dx + dy*dy + dz*dz)),
               -(dy*(-2*dx*dx + dy*dy + dz*dz)), -(dx*(dx*dx - 2*dy*dy + dz*dz)), (3*dx*dy*dz),
               -(dz*(-2*dx*dx + dy*dy + dz*dz)), (3*dx*dy*dz), -(dx*(dx*dx + dy*dy - 2*dz*dz)),
               -(dy*(-2*dx*dx + dy*dy + dz*dz)), -(dx*(dx*dx - 2*dy*dy + dz*dz)), (3*dx*dy*dz),
               -(dx*(dx*dx - 2*dy*dy + dz*dz)), -(3*dy*(dx*dx + dz*dz)), -(dz*(dx*dx - 2*dy*dy + dz*dz)),
               (3*dx*dy*dz), -(dz*(dx*dx - 2*dy*dy + dz*dz)), -(dy*(dx*dx + dy*dy - 2*dz*dz)),
               -(dz*(-2*dx*dx + dy*dy + dz*dz)), (3*dx*dy*dz), -(dx*(dx*dx + dy*dy - 2*dz*dz)),
               (3*dx*dy*dz), -(dz*(dx*dx - 2*dy*dy + dz*dz)), -(dy*(dx*dx + dy*dy - 2*dz*dz)),
               -(dx*(dx*dx + dy*dy - 2*dz*dz)), -(dy*(dx*dx + dy*dy - 2*dz*dz)),-(3*dz*(dx*dx + dy*dy));
    d2dfdp2 = d2dfdp2/(dnorm*dnorm*dnorm*dnorm*dnorm);
    MatrixXd Rstack(3,9);
    Rstack << (x.R*e).transpose(), MatrixXd::Zero(1,6),
              MatrixXd::Zero(1,3), (x.R*e).transpose(), MatrixXd::Zero(1,3),
              MatrixXd::Zero(1,6), (x.R*e).transpose();
    if (Lx){
      pop_Lx(x,h,Lx);
      //Analytical Lx
      //Matrix<double, _nx, 1> LxTrue;
      /*LxTrue.setZero();
      Vector3d dLdR = gain*(ehat.transpose()*x.R.transpose()*df);
      Vector3d dLdp = gain*(ddfdp.transpose()*x.R*e);
      LxTrue.head(3) = dLdR;
      LxTrue.segment(3,3) = dLdp;*/
      //Finite
      //Matrix<double, _nx, 1> LxF;
      //finite_Lx(x,h,&LxF);
      //if (!useFinite) {
      //  (*Lx) = LxTrue;
      //} else {
      //  (*Lx) = LxF;
      //}
      //cout << "LxA: " << endl << LxTrue << endl;
      //cout << "LxF: " << endl << LxF << endl;
      //cout << "Lx: " << endl << (*Lx) << endl;
    }
    if (Lxx){
      /*
      //Analytical
      Matrix<double, _nx, _nx> LxxTrue;
      LxxTrue.setZero();
      Matrix3d RTdfhat;
      SO3::Instance().hat(RTdfhat,x.R.transpose()*df);
      Matrix3d d2LdR2 = -gain*ehat*RTdfhat;
      Matrix3d d2LdRdp = -gain*(ehat.transpose()*x.R.transpose()*ddfdp);
      Matrix3d d2LdpdR = -gain*(ddfdp.transpose()*x.R*ehat);
      Matrix3d d2Ldp2 = -gain*(Rstack*d2dfdp2);
      LxxTrue.topLeftCorner(3,3) = d2LdR2;
      LxxTrue.block(0,3,3,3) = d2LdRdp;
      LxxTrue.block(3,0,3,3) = d2LdpdR;
      LxxTrue.block(3,3,3,3) = d2Ldp2;
      //Finite
      Matrix<double, _nx, _nx> LxxF;
      LxxF.setZero();
      double eps = 1e-10;
      Vectornd dx;
      for (int i = 0; i < Cost<T, _nx, _nu, _np>::sys.X.n; ++i ) {
        //Set dx to eps in the the i-th dim
        dx.setZero();
        dx[i] = eps;
        //Calculate using dx
        T xtemp;
        Cost<T, _nx, _nu, _np>::sys.X.Retract(xtemp,x,dx);//xtemp = xa+dx
        Cost<T, _nx, _nu, _np>::sys.Rec(xtemp,h);//Reconstruct
        Matrix<double, _nx, 1> Lpos;
        pop_Lx(xtemp,h, &Lpos);
        Lpos = Lpos - (*Lx);
        //Calculate using -dx
        Cost<T, _nx, _nu, _np>::sys.X.Retract(xtemp,x,-dx);//xtemp = xa+dx
        Cost<T, _nx, _nu, _np>::sys.Rec(xtemp,h);//Reconstruct
        Matrix<double, _nx, 1> Lneg;
        pop_Lx(xtemp,h, &Lneg);
        Lneg = Lneg - (*Lx);
        //Set Lx[i]
        LxxF.col(i) = (Lpos - Lneg)/(2*eps);
      }
      if(!useFinite) {
        (*Lxx) = LxxTrue;
      } else {
        (*Lxx) = LxxF;
      }
      cout << "LxxA: " << endl << LxxTrue << endl;
      cout << "LxxF: " << endl << LxxF << endl;
      cout << "Lxx: " << endl << (*Lxx) << endl;*/
      Lxx->setIdentity();
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

    //cout << "Yaw L Returned: " << c << endl;
    return c;
  }
}

#endif
