#ifndef GCOP_CONSTRAINTCOST_H
#define GCOP_CONSTRAINTCOST_H

#include <limits>
#include "lscost.h"
#include "constraint.h"
#include <iostream>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  template <typename T, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic,
    int _ng = Dynamic> class ConstraintCost : public LsCost<T, _nx, _nu, _np, _ng> {
    public:

  typedef Matrix<double, _ng, 1> Vectorgd;
  typedef Matrix<double, _ng, _nx> Matrixgnd;
  typedef Matrix<double, _ng, _nu> Matrixgcd;
  typedef Matrix<double, _ng, _np> Matrixgmd;

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
     * Linear-quadratic cost on a general manifold. Use this constructor for dynamic-size
     * control problem, i.e. ConstraintCost<T>(X, U, tf, xf, ...)
     * @param sys system
     * @param tf final time
     * @param con constraint
     */
    ConstraintCost(System<T, _nx, _nu, _np> &sys, double tf, Constraint<T, _nx, _nu, _np, _ng> &con);
    
    virtual double L(double t, const T& x, const Vectorcd& u, double h,
                     const Vectormd *p = 0,
                     Vectornd *Lx = 0, Matrixnd* Lxx = 0,
                     Vectorcd *Lu = 0, Matrixcd* Luu = 0,
                     Matrixncd *Lxu = 0,
                     Vectormd *Lp = 0, Matrixmd *Lpp = 0,
                     Matrixmnd *Lpx = 0);

    bool Res(Vectorgd &g, 
             double t, const T &x, const Vectorcd &u, double h,
             const Vectormd *p = 0,
             Matrixgnd *dgdx = 0, Matrixgcd *dgdu = 0,
             Matrixgmd *dgdp = 0);                 

    Constraint<T, _nx, _nu, _np, _ng> &con;

    double b;    ///< penalty coefficient (default is 1)

    protected:

    Vectorgd gp;      ///< active residuals
    Matrixgnd dgdx;   ///< gradients

    };
  
  template <typename T, int _nx, int _nu, int _np, int _ng>
    ConstraintCost<T, _nx, _nu, _np, _ng>::ConstraintCost(System<T, _nx, _nu, _np> &sys, 
                                                          double tf, Constraint<T, _nx, _nu, _np, _ng> &con) : 
    LsCost<T, _nx, _nu, _np, _ng>(sys, tf, con.ng), con(con), b(1) {

    if (_ng == Dynamic) {
      gp.resize(con.ng);
      this->ng = con.ng;
    }    
    else 
    {
      this->ng = _ng;
    }

    if (_nx == Dynamic || _ng == Dynamic) {
      dgdx.resize(con.ng, sys.X.n);
    }
  }

  template <typename T, int _nx, int _nu, int _np, int _ng> 
    bool  ConstraintCost<T, _nx, _nu, _np, _ng>::Res(Vectorgd &g, 
             double t, const T &x, const Vectorcd &u, double h,
             const Vectormd *p,
             Matrixgnd *dgdx, Matrixgcd *dgdu,
             Matrixgmd *dgdp)
  {
    this->con(this->g, t, x, u, p, dgdx);
    for (int i = 0; i < con.ng; ++i) {
      if (this->g[i] < 0) { // if constraint is satisfied then null it
        this->g[i] = 0;
        if (dgdx)
          dgdx->row(i).setZero();
      }      
    }
    g = sqrt(b)*this->g;
    return true;
  }
  
  template <typename T, int _nx, int _nu, int _np, int _ng> 
    double ConstraintCost<T, _nx, _nu, _np, _ng>::L(double t, const T &x, const Matrix<double, _nu, 1> &u,
                                                    double h,
                                                    const Matrix<double, _np, 1> *p,
                                                    Matrix<double, _nx, 1> *Lx, Matrix<double, _nx, _nx> *Lxx,
                                                    Matrix<double, _nu, 1> *Lu, Matrix<double, _nu, _nu> *Luu,
                                                    Matrix<double, _nx, _nu> *Lxu,
                                                    Matrix<double, _np, 1> *Lp, Matrix<double, _np, _np> *Lpp,
                                                    Matrix<double, _np, _nx> *Lpx) {

    double q = 234234023411230; 
    dgdx(0,0) = q;   // random number

    // only consider state constraints for now
    this->con(this->g, t, x, u, p, Lx ? &dgdx : 0);
    
    T xb;
    Vectornd dx;
    Vectorgd gp;
    Vectorgd gm;
    double eps = 1e-3;

    // if no jacobians were provided use finite differences
    if (Lx && fabs(dgdx(0,0) - q) < 1e-10) {
      
      for (int i = 0; i < this->sys.X.n; ++i) {
        dx.setZero();
        dx[i] = eps;
            
        // xb = x + dx
        this->sys.X.Retract(xb, x, dx);
        
        // reconstruct state using previous time-step
        this->sys.Rec(xb, h);
            
        // gp = con(xb)
        this->con(gp, t, xb, u, p);
        
        // xb = x - dx
        this->sys.X.Retract(xb, x, -dx);
            
        // reconstruct state using previous time-step
        this->sys.Rec(xb, h);

        // gm = con(xb)
        this->con(gm, t, xb, u, p);

        dgdx.col(i) = (gp - gm)/(2*eps);
      }
    }

    for (int i = 0; i < con.ng; ++i) {
      if (this->g[i] < 0) { // if constraint is satisfied then null it
        this->g[i] = 0;
        if (Lx)
          dgdx.row(i).setZero();
      }      
    }

    if (Lx)
      *Lx = b*dgdx.transpose()*this->g;
    if (Lxx) {
      *Lxx = b*dgdx.transpose()*dgdx;  // use a GN approximation to the Hessian
      
      /*
      if (this->g[0] > 0) {
        Vector3d v = dgdx.segment(3,3);
        Matrix3d M= Matrix3d::Identity() - v*v.transpose();
        Lxx->block(3,3,3,3) += -b*M;
      }
      */
    }
    /*    if (Lxx) {
      Lxx->Identity();
      *Lxx = .01**Lxx;
      }
    */

    return b/2*this->g.squaredNorm();
    //return b/2*this->g.squaredNorm();
  }  
}


#endif
