#ifndef GCOP_LSCOST_H
#define GCOP_LSCOST_H

#include <Eigen/Dense>
#include <iostream>
#include "system.h"
#include "cost.h"

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
  template <typename T = VectorXd, 
    int _nx = Dynamic,
    int _nu = Dynamic,
    int _np = Dynamic,
    int _ng = Dynamic> class LsCost : public Cost<T, _nx, _nu, _np> {
    
  public:
  
  typedef Matrix<double, _ng, 1> Vectorgd;
  typedef Matrix<double, _ng, _nx> Matrixgxd;
  typedef Matrix<double, _ng, _nu> Matrixgud;
  typedef Matrix<double, _ng, _np> Matrixgpd;
  
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
   * _nureate a _nuost interfacex
   * @param X the state manifold which is used to perform addition/subtraction of states
   * @param tf time horizon: when the cost function L is called it will internally check whether its argument t is equation to tf and return the terminal cost
   * @param ng number of residuals per time-step
   */
  LsCost(System<T, _nx, _nu, _np> &sys, double tf, int ng = 0);
  
  /**
   * Residual cost
   * @param g residual
   * @param t time
   * @param x state
   * @param u control
   * @param h time-step
   * @param p parameter (optional)
   * @param dgdx derivative wrt x
   * @param dgdu derivative wrt u
   * @param dgdp derivative wrt p
   */
  virtual bool Res(Vectorgd &g, 
                   double t, const T &x, const Vectorcd &u, double h,
                   const Vectormd *p = 0,
                   Matrixgxd *dgdx = 0, Matrixgud *dgdu = 0,
                   Matrixgpd *dgdp = 0);    

  int ng;       ///< number of residuals per time-step
  Vectorgd g;   ///< residuals per time-step
  };
  
  
  template <typename T, int _nx, int _nu, int _np, int _ng> 
    LsCost<T, _nx, _nu, _np, _ng>::LsCost(System<T, _nx, _nu, _np> &sys, 
                                          double tf, int ng) : 
    Cost<T, _nx, _nu, _np>(sys, tf), ng(_ng != Dynamic ? _ng : ng)
    {
      assert(ng > 0);
      if (_ng == Dynamic)
        g.resize(ng);
    }
  
  template <typename T, int _nx, int _nu, int _np, int _ng> 
    bool LsCost<T, _nx, _nu, _np, _ng>::Res(Vectorgd &g, 
                                            double t, const T &x, const Vectorcd &u, double h,
                                            const Vectormd *p,
                                            Matrixgxd *dgdx, Matrixgud *dgdu,
                                            Matrixgpd *dgdp) {
    cout << "[W] LsCost:Res: unimplemented! Subclasses should override." << endl;
    return false;
  }
}

#endif

