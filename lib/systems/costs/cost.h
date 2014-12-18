#ifndef GCOP_COST_H
#define GCOP_COST_H

#include <Eigen/Dense>
#include <iostream>
#include "system.h"

namespace gcop {

  using namespace Eigen;
  using namespace std;
  
  /**
   * Cost interface for optimal control on manifolds. 
   * Defines a cost function and means to compute difference between
   * two states on a manifold. 
   *
   * Subclasses should provide implementation of either a regular const function L
   * or a parameter-dependent cost function Lp 
   * (e.g. for sys id / adaptive control / optimal design problems)
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <typename T, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic> class Cost {
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
  Cost(System<T, _nx, _nu, _np> &sys, double tf);

  /**
   * Create a Sensor Based Cost interface
   * @param sys the system for which cost function is defined. Provides the system manifold X
   * @param  
   */
  //Cost(System<T, _nx, _nu, _np> &sys, Sensor<double tf);
  
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
  
  System<T, _nx, _nu, _np> &sys;  ///< system  
  double tf;               ///< final time of trajectory
};


  template <typename T, int _nx, int _nu, int _np> 
    Cost<T, _nx, _nu, _np>::Cost(System<T, _nx, _nu, _np> &sys, 
                                 double tf) : sys(sys), tf(tf) {
  }  
  
  template <typename T, int _nx, int _nu, int _np> 
    double Cost<T, _nx, _nu, _np>::L(double t, const T& x, const Vectorcd& u, double h,
                                     const Vectormd *p,
                                     Matrix<double, _nx, 1> *Lx, Matrix<double, _nx, _nx>* Lxx,
                                     Matrix<double, _nu, 1> *Lu, Matrix<double, _nu, _nu>* Luu,
                                     Matrix<double, _nx, _nu> *Lxu,
                                     Vectormd *Lp, Matrixmd *Lpp,
                                     Matrixmnd *Lpx) {
    cout << "[W] Cost:L: unimplemented!" << endl;
    return 0;
  }
}

#endif
