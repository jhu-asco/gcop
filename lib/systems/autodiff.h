#ifndef GCOP_AUTODIFF_H
#define GCOP_AUTODIFF_H

#include <Eigen/Dense>
#include <vector>
#include <assert.h>
#include "system.h"
#include <iostream>

#define NUMBER_DIRECTIONS 16
#include <unsupported/Eigen/AdolcForward>

int adtl::ADOLC_numDir;

namespace gcop {
  
  using namespace Eigen;
  
  template <typename T, int _n = Dynamic, int _c = Dynamic> 
    class Autodiff {
    
  public:
  
  typedef Matrix<double, _n, 1> Vectornd;
  typedef Matrix<double, _c, 1> Vectorcd;
  typedef Matrix<double, _n, _n> Matrixnd;
  typedef Matrix<double, _n, _c> Matrixncd;
  typedef Matrix<double, _c, _n> Matrixcnd;
  typedef Matrix<double, _c, _c> Matrixcd;
  
  typedef Matrix<double, Dynamic, 1> Vectormd;
  typedef Matrix<double, Dynamic, Dynamic> Matrixmd;
  typedef Matrix<double, _n, Dynamic> Matrixnmd;
  typedef Matrix<double, Dynamic, _n> Matrixmnd;
  

  Autodiff(System<T, Vectorcd, _n, _c> &sys);
  
  double DF(Vectornd &v, double t, const T& xa, 
            const Vectorcd &u, double h,
            Matrixnd *A, Matrixncd *B);
  
  double DF(Vectornd &v, double t, const T& xa,
            const Vectorcd &u, double h,
            const Vectormd &p,
            Matrixnd *A, Matrixncd *B, Matrixnmd *C);


  struct Fdx {
  Fdx(System<T, Vectorcd, _n, _c> &sys) : sys(sys) {};
    
    typedef Vectornd InputType;
    typedef Vectornd ValueType;
    typedef Matrixnd JacobianType;
    
    //    Fdx(double t, const T &x, const Vectorcd& u
    //    int inputs() const { return n; }
    //    int values() const { return n; }
    
    void operator() (const Vectornd& e, Vectornd* v) const {
      sys.X.Retract((T&)xe, *x, e);
      sys.F(*v, t, xe, *u, h);
    }
    
    System<T, Vectorcd, _n, _c> &sys;
    double t;
    T xe;
    const T *x;
    const Vectorcd* u;
    double h;
  };
  
  struct Fu {      
  Fu(System<T, Vectorcd, _n, _c> &sys) : sys(sys) {};

    typedef Vectorcd InputType;
    typedef Vectornd ValueType;
    typedef Matrixncd JacobianType;

    //    int inputs() const { return c; }
    //    int values() const { return n; }
    
    void operator() (const Vectorcd& u, Vectornd* v) const {
      sys.F(*v, t, *x, u, h);
    }

    System<T, Vectorcd, _n, _c> &sys;
    double t;
    const T* x;
    double h;
  };
  
  
  Fdx fdx; 
  Fu fu;
  
  };
  
  
  
  template <typename T, int _n, int _c> 
    Autodiff<T, _n, _c>::Autodiff(System<T, Matrix<double, _c, 1>, _n, _c> &sys) : 
    fdx(sys), fu(sys) {
    adtl::ADOLC_numDir = NUMBER_DIRECTIONS;
  }

  
  template <typename T, int _n, int _c> 
    double Autodiff<T, _n, _c>::DF(Vectornd &v, double t, const T &x, 
                                   const Matrix<double, _c ,1> &u, double h,
                                   Matrix<double, _n, _n> *A, Matrix<double, _n, _c> *B) {
    

    assert(A);
    assert(B);

    Matrixnd Dfdx;
    Matrixnd Ad;
    Matrixnd D;
    Vectornd e;

    
    fdx.t = t; 
    fdx.x = &x;
    fdx.xe = x;
    fdx.u = &u;
    fdx.h = h;
    AdolcForwardJacobian<System::Fdx> autoj(fdx);
    if (_n == Dynamic) {
      Dfdx.resize(n,n);
      Ad.resize(n,n);
      D.resize(n,n);
      e.resize(n);
    }
    
    e.setZero();
    
    autoj(e, &v, &Dfdx);
    
    X.Adtau(Ad, -v);
    X.dtau(D, -v);
    
    *A = Ad + D*Dfdx;
    
    // B matrix
    Matrixncd Dfu;
    if (_n == Dynamic || _c == Dynamic)
      Dfu.resize(n, c);
    
    fu.t = t; 
    fu.x = &xa;
    fu.h = h;
    AdolcForwardJacobian<Fdx> autoju(fu);
     
    autoju(u, &v, &Dfu);
    
    *B = D*Dfu;
    
    return 0;    
  }
  
  template <typename T, int _n, int _c> 
    double Autodiff<T, _n, _c>::F(Vectornd &v, double t, const T& xa,
                                  const Matrix<double, _c, 1> &u, double h,
                                  const VectorXd &p,
                                  Matrix<double, _n, _n> *A, Matrix<double, _n, _c> *B, 
                                  Matrix<double, _n, Dynamic> *C) {
    
  }
}

#endif
