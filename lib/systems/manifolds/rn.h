#ifndef GCOP_RN_H
#define GCOP_RN_H

#include "manifold.h"
#include <cfloat>

namespace gcop {
  
  using namespace Eigen;
  
  /**
   * n-dimension vector space, i.e. \f$ X = R^n \f$
   */
  template <int _n = Dynamic> class Rn :
    public Manifold< Matrix<double, _n, 1>, _n> {
    
    typedef Matrix<double, _n, 1> Vectornd;
    typedef Matrix<double, _n, _n> Matrixnd;
    
  public:

    Rn(int n = _n);
    
    /**
     * Singleton to be used only for fixed-size Rn, i.e. Rn<n>::Instance()
     * @return an instance of Rn<n>
     */
    static Rn<_n>& Instance();
    
    void Lift(Vectornd &v,
              const Vectornd &xa,
              const Vectornd &xb);
    
    void Retract(Vectornd &xb,
                 const Vectornd &xa,
                 const Vectornd &v);

    void dtau(Matrixnd &M, const Vectornd &v);
    
    void Adtau(Matrixnd &M, const Vectornd &v);

    static const char BND_BOX = 0;
    static const char BND_ELLIPSOID = 1;    

    char bndType; ///< type of bound
    Vectornd w;   ///< ellipsoidal bound weights

    /**
     * Find a vector v = v0 + d for given v0 and d  
     * such that v satisfies the bounds; this is accomplished by 
     * modifying d appropriately and returning the updated v and d.
     * Default implementation is assuming the bound is a box defined by
     * the vectors lb and ub.
     * @param v vector satisfying bounds
     * @param d direction
     * @param v0 starting vector
     * @return true on success, false if the arguments are not compatible
     */
    virtual bool Bound(Vectornd &v, Vectornd &d, const Vectornd &v0);

  };

  template <int _n> 
    Rn<_n>::Rn(int n) : Manifold<Matrix<double, _n, 1>, _n>(n), bndType(BND_BOX) {
    if (_n == Dynamic) {
      this->lb.resize(n);
      this->ub.resize(n);
      w.resize(n);
    }
    this->lb.setConstant(-1e16);
    this->ub.setConstant(1e16);      
    w.setConstant(1);
  }
  
  template <int _n> 
    Rn<_n>& Rn<_n>::Instance() {
    static Rn<_n> instance;
    return instance;
  }

  template <int _n>
    void Rn<_n>::Lift(Matrix<double, _n, 1> &v, 
                     const Matrix<double, _n, 1> &xa,
                     const Matrix<double, _n, 1> &xb) {
    v = xb - xa;
  }

  template <int _n> 
    void Rn<_n>::Retract(Matrix<double, _n, 1> &xb, 
                         const Matrix<double, _n, 1> &xa,
                         const Matrix<double, _n, 1> &v) {
    xb = xa + v;
  }

  template <int _n> 
    void Rn<_n>::dtau(Matrix<double, _n, _n> &M, 
                      const Matrix<double, _n, 1> &v) {
    M.setIdentity();
  }

  template <int _n> 
    void Rn<_n>::Adtau(Matrix<double, _n, _n> &M, 
                       const Matrix<double, _n, 1> &v) {
    M.setIdentity();
  }

  template <int _n> 
    bool Rn<_n>::Bound(Vectornd &v, Vectornd &d, const Vectornd &v0) {

    if (bndType == BND_BOX) {
      for (int j = 0; j < v.size(); ++j) {
        if (v[j] < this->lb[j]) {
          v[j] = this->lb[j];
          d[j] = this->lb[j] - v0[j];
        } else
          if (v[j] > this->ub[j]) {
            v[j] = this->ub[j];
            d[j] = this->ub[j] - v0[j];
        }
      } 
      return true;
    }

    if (bndType == BND_ELLIPSOID) {
      
      // if constraint satisfied
      if ((v.cwiseProduct(w)).dot(v) <= 1)
        return false;

      // if constraint cannot be satisfied
      if ((v0.cwiseProduct(w)).dot(v0) >= 1) {
        v = v0;
        d.setZero();
        return false;        
      }

      double a = (d.cwiseProduct(w)).dot(d);
      double b = 2*(d.cwiseProduct(w)).dot(v0);
      double c = (v0.cwiseProduct(w)).dot(v0) - 1; // should be negative
      assert(c <= 1e-16);

      double s = (- b + sqrt(b*b - 4*a*c))/(2*a);
      
      //      assert(s >= -1e-16);
      
      // p
      for (int i = 0; i < d.size(); ++i) {
        if (w[i] > 1e-16) {
          d[i] = s*d[i];   
          v[i] = v0[i] + d[i];
        }        
      }
      return true;
    }
    
  }

}

#endif
