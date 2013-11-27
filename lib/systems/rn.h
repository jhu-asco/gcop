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
  };

  template <int _n> 
    Rn<_n>::Rn(int n) : Manifold<Matrix<double, _n, 1>, _n>(n) {
    if (_n == Dynamic) {
      this->lb.resize(n);
      this->ub.resize(n);
    }
    this->lb.setConstant(-1e16);
    this->ub.setConstant(1e16);      
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

}

#endif
