#ifndef GCOP_MULTILSCOST_H
#define GCOP_MULTILSCOST_H

#include <Eigen/Dense>
#include <iostream>
#include "lscost.h"

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
    int _np = Dynamic> class MultiLsCost : public LsCost<T, _nx, _nu, _np> {
  public:

  typedef Matrix<double, Dynamic, 1> Vectorgd;
  typedef Matrix<double, Dynamic, _nx> Matrixgxd;
  typedef Matrix<double, Dynamic, _nu> Matrixgud;
  typedef Matrix<double, Dynamic, _np> Matrixgpd;
  
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
  MultiLsCost(System<T, _nx, _nu, _np> &sys, double tf);
  
  /**
   *  Add a cost
   * @param cost the cost
   */
  virtual void AddCost(LsCost<T, _nx, _nu, _np>* cost);
  
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
  
  Vectornd Lx;
  Matrixnd Lxx;
  Vectorcd Lu;
  Matrixcd Luu;  
  Matrixncd Lxu;
  Vectormd Lp;
  Matrixmd Lpp;
  Matrixmnd Lpx;

protected:

  vector<LsCost<T, _nx, _nu, _np>* > costs;
  };


  template <typename T, int _nx, int _nu, int _np> 
    MultiLsCost<T, _nx, _nu, _np>::MultiLsCost(System<T, _nx, _nu, _np> &sys, 
                                           double tf) : LsCost<T, _nx, _nu, _np>(sys, tf, 1) {
    this->ng = 0;
    if (_nx == Dynamic) {
      Lx.resize(sys.X.n);
      Lxx.resize(sys.X.n, sys.U.n);
      Lxu.resize(sys.X.n, sys.U.n);
      Lpx.resize(sys.P.n, sys.X.n);      
    }
    if (_nu == Dynamic) {
      Lu.resize(sys.U.n);
      Luu.resize(sys.U.n, sys.U.n);      
    }
    if (_np == Dynamic) {
      Lp.resize(sys.P.n);
      Lpp.resize(sys.P.n, sys.P.n);      
    }
  }  

  template <typename T, int _nx, int _nu, int _np> 
    void MultiLsCost<T, _nx, _nu, _np>::AddCost(LsCost<T, _nx, _nu, _np>* cost)
  {
    this->ng += cost->ng;
    this->g.resize(this->ng);
    costs.push_back(cost);
  }

  template <typename T, int _nx, int _nu, int _np> 
    bool MultiLsCost<T, _nx, _nu, _np>::Res(Vectorgd &g, double t, const T& x, const Vectorcd& u, 
                                          double h,
                                          const Vectormd *p,
                                          Matrixgxd *dgdx, Matrixgud *dgdu, Matrixgpd *dgdp)
  {
    g = Vectorgd::Zero(this->ng);
    if (dgdx)
      dgdx->setZero();
    if (dgdu)
      dgdu->setZero();
    if (dgdp)
      dgdp->setZero();
    int start_idx = 0;

    for (int i = 0; i < costs.size(); ++i) {
      //Matrix<double, Dynamic, _nx> dgdx_i(costs[i]->ng, _nx);
      //Matrix<double, Dynamic, _nu> dgdu_i(costs[i]->ng, _nu);
      //Matrix<double, Dynamic, _np> dgdp_i(costs[i]->ng, _np);
      VectorXd g_seg;
      costs[i]->Res(g_seg, t, x, u, h, p, 0, 0, 0);
      //                    (dgdx ? &dgdx_i : 0),
      //                    (dgdu ? &dgdu_i : 0), 
      //                    dgdp ? &dgdp_i : 0);
      g.segment(start_idx, costs[i]->ng) = g_seg;
      /*
      if (dgdx)
        dgdx->block(start_idx, 0, costs[i]->ng, _nx) = dgdx_i;
      if (dgdu)
        dgdu->block(start_idx, 0, costs[i]->ng, _nu) = dgdu_i;
      if (dgdp)
        dgdp->block(start_idx, 0, costs[i]->ng, _np) = dgdp_i;
      */
      start_idx += costs[i]->ng;
    }
    return true;
  }
  
  template <typename T, int _nx, int _nu, int _np> 
    double MultiLsCost<T, _nx, _nu, _np>::L(double t, const T& x, const Vectorcd& u, double h,
                                          const Vectormd *p,
                                          Matrix<double, _nx, 1> *Lx, Matrix<double, _nx, _nx>* Lxx,
                                          Matrix<double, _nu, 1> *Lu, Matrix<double, _nu, _nu>* Luu,
                                          Matrix<double, _nx, _nu> *Lxu,
                                          Vectormd *Lp, Matrixmd *Lpp,
                                          Matrixmnd *Lpx) {
    double c = 0;
    if (Lx)
      Lx->setZero();
    if (Lxx)
      Lxx->setZero();
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

    for (int i = 0; i < costs.size(); ++i) {
      c = c + costs[i]->L(t, x, u, h, p,
                          (Lx ? &this->Lx : 0),
                          (Lxx ? &this->Lxx : 0), 
                          Lu ? &this->Lu : 0, 
                          Luu ? &this->Luu : 0,
                          Lxu ? &this->Lxu : 0, 
                          Lp ? &this->Lp : 0, 
                          Lpp ? &this->Lpp : 0, 
                          Lpx ? &this->Lpx : 0);
      if (Lx)
        *Lx = *Lx + this->Lx;
      if (Lxx)
        *Lxx = *Lxx + this->Lxx;
      if (Lu)
        *Lu = *Lu + this->Lu;
      if (Luu)
        *Luu = *Luu + this->Luu;
      if (Lxu)
        *Lxu = *Lxu + this->Lxu;
      if (Lp)
        *Lp = *Lp + this->Lp;
      if (Lpp)
        *Lpp = *Lpp + this->Lpp;
      if (Lpx)
        *Lpx = *Lpx + this->Lpx;
    }
    return c;
  }
}

#endif
