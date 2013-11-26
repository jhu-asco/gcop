#ifndef GCOP_FUNCTION_H
#define GCOP_FUNCTION_H

#include <Eigen/Dense>
#include "manifold.h"
#include "rn.h"
#include <assert.h>

namespace gcop {
  
  using namespace Eigen;
  /**
   * m-dimensional vector-valued function on an n-dimensional manifold X
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <typename T, int _n = Dynamic, int _m = Dynamic> class Function {
  public:
    typedef Matrix<double, _n, 1> Vectornd;
    typedef Matrix<double, _m, 1> Vectormd;
    typedef Matrix<double, _m, _n> Matrixmnd;
    
    Function(Manifold<T, _n> &X, double eps = 1e-12);

    virtual void F(Vectormd &f, const T &x) = 0;

    void DF(Matrixmnd &Df, const T &x);    

    Manifold<T, _n> &X;
    double eps;   ///< ste-size of finite differences
  };

  
  template <typename T, int _n, int _m>
    Function<T, _n, _m>::Function(Manifold<T, _n> &X, double eps) : 
    X(X),
    eps(eps) {    
    }
  

  template <typename T, int _n, int _m>
    void Function<T, _n, _m>::DF(Matrix<double, _m, _n> &Df, 
                                 const T &x) {

    typedef Matrix<double, _n, 1> Vectornd;
    typedef Matrix<double, _m, 1> Vectormd;
    typedef Matrix<double, _m, _n> Matrixmnd;
    
    int m = Df.rows();
    int n = Df.cols();

    assert(n > 0);
    assert(m > 0);

    T xl(x);
    T xr(x);
    Vectormd fl;
    Vectormd fr;
    Vectornd e;

    if (_m == Dynamic) {
      fl.resize(m);
      fr.resize(m);
    } else {
      assert(_m == m);
    }

    if (_n == Dynamic) {
      e.resize(n);
    } else {
      assert(_n == n);
    }
    
    for (int i = 0; i < n; ++i) {      
      e.setZero();
      
      // left
      e[i] = -eps;
      X.Retract(xl, x, e);
      F(fl, xl);

      // right
      e[i] = eps;
      X.Retract(xr, x, e);
      F(fr, xr);

      // central difference
      Df.block(0, i, m, 1) = (fr - fl)/(2*eps);                           
    }
  }
  
}

#endif

