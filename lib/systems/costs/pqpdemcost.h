#ifndef GCOP_PQPDEMCOST_H
#define GCOP_PQPDEMCOST_H

#include <limits>
#include "constraintcost.h"
#include "pqpdem.h"
#include <iostream>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  template <typename T, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic> class PqpDemCost : public ConstraintCost<T, _nx, _nu, _np, 1> {
    public:

  typedef Matrix<double, 1, 1> Vectorgd;
  typedef Matrix<double, 1, _nx> Matrixgxd;
  typedef Matrix<double, 1, _nu> Matrixgud;
  typedef Matrix<double, 1, _np> Matrixgpd;

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
     * control problem, i.e. PqpdemCost<T>(X, U, tf, xf, ...)
     * @param sys system
     * @param tf final time
     * @param con pqpdem
     */
    PqpDemCost(System<T, _nx, _nu, _np> &sys, tf, PqpDem<T, _nx, _nu, _np, 1> &con);   

    };
  
  template <typename T, int _nx, int _nu, int _np>
    PqpDemCost<T, _nx, _nu, _np>::PqpdemCost(System<T, _nx, _nu, _np> &sys, 
                                             tf, PqpDem<T, _nx, _nu, _np> &con) : 
    ConstraintCost(sys, tf, con)
    LsCost<T, _nx, _nu, _np, 1>(sys, tf, con.ng), con(con), b(1) {

    if (1 == Dynamic) {
      gp.resize(con.gn);
    }    

    if (_nx == Dynamic || 1 == Dynamic) {
      dgdx.resize(con.gn, sys.X.n);
    }
  }
  
  template <typename T, int _nx, int _nu, int _np, int 1> 
    double PqpdemCost<T, _nx, _nu, _np, 1y>::L(double t, const T &x, const Matrix<double, _nu, 1> &u,
                                                    double h,
                                                    const Matrix<double, _np, 1> *p,
                                                    Matrix<double, _nx, 1> *Lx, Matrix<double, _nx, _nx> *Lxx,
                                                    Matrix<double, _nu, 1> *Lu, Matrix<double, _nu, _nu> *Luu,
                                                    Matrix<double, _nx, _nu> *Lxu,
                                                    Matrix<double, _np, 1> *Lp, Matrix<double, _np, _np> *Lpp,
                                                    Matrix<double, _np, _nx> *Lpx) {

    // only consider state pqpdems for now
    this->con(this->g, t, x, u, p, Lx ? dgdx : 0);

    for (int i=0; i<con.ng; ++i) {
      if (this->g[i] < 0) { // if pqpdemy is satisfied then null it
        this->g[i] = 0;
        if (Lx)
          dgdx.row(i).setZero();
      }      
    }
    if (Lx)
      *Lx = b*dgdx.transpose()*this->g;
    if (Lxx)
      *Lxx = b*dgdx.transpose()*dgdx;

    return b/2*this->g.squaredNorm();
  }  
}


#endif
