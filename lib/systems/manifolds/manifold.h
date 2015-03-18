#ifndef GCOP_MANIFOLD_H
#define GCOP_MANIFOLD_H

#include <Eigen/Dense>
#include <assert.h>
#include <iostream>

namespace gcop {
  
  using namespace Eigen;
  /**
   * Homogeneous manifold
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <typename T, int _n = Dynamic> class Manifold {
  public:
  typedef Matrix<double, _n, 1> Vectornd;
  typedef Matrix<double, _n, _n> Matrixnd;

  /**
   * Create a manifold of dimension n. For fixed-size use default constructor Manifold<T,n>()
   * For dynamic-size use Manifold<T>(n).
   *
   * @param n manifold dimension
   */
  Manifold(int n = _n);
  

  /**
   * Compute a lie algebra element corresponding to a group action between 
   * two given states on a manifold. This operation geometrically means
   * that the curve from xa to xb is "lifted" to the tangent space of the manifold
   * to a vector represented by the lie algebra element.
   *
   * The element can also be regarded as a difference between the 
   * two states xa and xb. On a vector space
   * this is v = xb - xa. On homogeneous manifolds v is a Lie algebra 
   * element corresponding to an action taking xa to xb. 
   *
   * @param v Lie algebra element (regarded as difference between xb and xa)
   * @param xa starting state
   * @param xb ending state
   * 
   */
  virtual void Lift(Vectornd &v, const T &xa, const T &xb) = 0;
  
  
  /**
   * Retracts a tangent vector on a manifold at state xa to produce a new state xb.
   * The vector is represented using a Lie algebra element.
   *
   * The operation can also be regarded as a way to perform addition on manifolds, i.e. 
   * to increment a given state xa by v to get a new state xb. On a vector space
   * this is xb = xa + v. On homogeneous manifolds v is a Lie algebra 
   * element corresponding to an action taking xa to xb.
   *
   * @param xb ending state
   * @param xa starting state
   * @param v Lie algebra element (regarded as difference between xb and xa)
   */
  virtual void Retract(T &xb, const T &xa, const Vectornd &v) = 0;


  /**
   * Right-trivialized derivative of the retraction map dtau(v)
   *
   * @param M matrix operator
   * @param v Lie algebra element
   */
  virtual void dtau(Matrixnd &M, const Vectornd &v) { M.setIdentity(); };

  /**
   * Adjoint map of the retraction, i.e. Ad(\tau(v))
   *
   * @param M resulting matrix operator
   * @param v Lie algebra element
   */
  virtual void Adtau(Matrixnd &M, const Vectornd &v) {M.setIdentity(); };
   
  int n;  ///< dimension
  
  bool bnd;    ///< is the space bounded (false by default). Bounds are defined below: these "box" bounds make sense when there is a natural ordering on the manifold, e.g. of the natural numbers)
  T lb;        ///< lower bound (-inf by default)
  T ub;        ///< upper bound (inf by default)

  };

  template <typename T, int _n>
    Manifold<T, _n>::Manifold(int n) : n(_n != Dynamic ? _n : n), bnd(false) {
  }
}

#endif
