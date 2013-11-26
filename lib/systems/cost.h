#ifndef GCOP_COST_H
#define GCOP_COST_H

#include <Eigen/Dense>
#include "manifold.h"

namespace gcop {

  using namespace Eigen;
  
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
  template <typename T, typename Tu, int _n = Dynamic, int _c = Dynamic> class Cost {
  public:
  
  typedef Matrix<double, _n, 1> Vectornd;
  typedef Matrix<double, _c, 1> Vectorcd;
  typedef Matrix<double, _n, _n> Matrixnd;
  typedef Matrix<double, _n, _c> Matrixnc;
  typedef Matrix<double, _c, _n> Matrixcn;
  typedef Matrix<double, _c, _c> Matrixcd;
  
  typedef Matrix<double, Dynamic, 1> Vectormd;
  typedef Matrix<double, Dynamic, Dynamic> Matrixmd;
  typedef Matrix<double, _n, Dynamic> Matrixnmd;
  typedef Matrix<double, Dynamic, _n> Matrixmnd;
  
  /**
   * _create a _cost interface
   * @param X the state manifold which is used to perform addition/subtraction of states
   * @param tf time horizon: when the cost function L is called it will internally check whether its argument t is equation to tf and return the terminal cost
   */
  Cost(Manifold<T, _n> &X, Manifold<Tu, _c> &U, double tf);
  
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
  virtual double L(double t, const T& x, const Tu& u, double h,
                   Vectornd *Lx = 0, Matrixnd* Lxx = 0,
                   Vectorcd *Lu = 0, Matrixcd* Luu = 0,
                   Matrixnc *Lxu = 0);
  
  
  /**
   * Parameter-dependent cost function L
   * @param t time
   * @param x state
   * @param u control
   * @param p parameter
   * @param Lx derivative wrt x
   * @param Lxx derivative wrt x,x
   * @param Lu derivative wrt u
   * @param Luu derivative wrt u,u
   * @param Lxu derivative wrt x,u     
   * @param Lp derivative wrt p
   * @param Lpp derivative wrt p,p
   * @param Lpx derivative wrt p,x     
   */
  virtual double Lp(double t, const T& x, const Tu& u,
                    const Vectormd &p,
                    Vectornd *Lx = 0, Matrixnd* Lxx = 0,
                    Vectorcd *Lu = 0, Matrixcd* Luu = 0,
                    Matrixnc *Lxu = 0,
                    Vectormd *Lp = 0, Matrixmd *Lpp = 0,
                    Matrixmnd *Lpx = 0);
  
  Manifold<T, _n> &X;   ///< the state manifold
  Manifold<Tu, _c> &U;  ///< the control manifold
  
  double tf;         ///< final time of trajectory
  };
  
  
  template <typename T, typename Tu, int _n, int _c> 
    Cost<T, Tu, _n, _c>::Cost(Manifold<T, _n> &X, Manifold<Tu, _c> &U, double tf) : X(X), U(U), tf(tf) {
  }  
  
  template <typename T, typename Tu, int _n, int _c> 
    double Cost<T, Tu, _n, _c>::L(double t, const T& x, const Tu& u, double h,
                                  Matrix<double, _n, 1> *Lx, Matrix<double, _n, _n>* Lxx,
                                  Matrix<double, _c, 1> *Lu, Matrix<double, _c, _c>* Luu,
                                  Matrix<double, _n, _c> *Lxu) {
    return 0;
  }
  
  template <typename T, typename Tu, int _n, int _c> 
    double Cost<T, Tu, _n, _c>::Lp(double t, const T &x, const Tu &u,
                                   const VectorXd &p,
                                   Matrix<double, _n, 1> *Lx, Matrix<double, _n, _n> *Lxx,
                                   Matrix<double, _c, 1> *Lu, Matrix<double, _c, _c> *Luu,
                                   Matrix<double, _n, _c> *Lxu,
                                   Matrix<double, Dynamic, 1> *Lp, Matrix<double, Dynamic, Dynamic> *Lpp,
                                   Matrix<double, Dynamic, _n> *Lpx) {
    return 0;
  }
}

#endif

