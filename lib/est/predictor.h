#ifndef GCOP_PREDICTOR_H
#define GCOP_PREDICTOR_H

#include "system.h"

namespace gcop {
  
  template <typename T = VectorXd, 
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic> class Predictor {
  public:

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;

  typedef Matrix<double, _nx, _nx> Matrixnd;
  typedef Matrix<double, _nx, _nu> Matrixncd;
  typedef Matrix<double, _nu, _nx> Matrixcnd;
  typedef Matrix<double, _nu, _nu> Matrixcd;
  
  typedef Matrix<double, _np, _np> Matrixmd;
  typedef Matrix<double, _nx, _np> Matrixnmd;
  typedef Matrix<double, _np, _nx> Matrixmnd;

  /**
   * Basic Predictor
   * @param sys system model
   */
  Predictor(System<T, _nx, _nu, _np>  &sys);
  
  virtual ~Predictor();
    
  /**
   * Prediction step. 
   * @param xb new belief state
   * @param t time
   * @param xa previous belief state
   * @param u control inputs
   * @param h time-step
   * @param p parameters (optional)
   * @param cov whether to update the covariance as well (true by default)
   * @return true if success
   */
  virtual bool Predict(T& xb, double t, const T &xa, 
                       const Vectorcd &u, double h, 
                       const Vectormd *p = 0, bool cov = true) = 0;
  
  
  System<T, _nx, _nu, _np>  &sys;             ///< system    

  };
  
  
  template <typename T, int _nx, int _nu, int _np> 
    Predictor<T, _nx, _nu, _np>::Predictor( System<T, _nx, _nu, _np>  &sys) : 
    sys(sys) {
  }  

  template <typename T, int _nx, int _nu, int _np> 
    Predictor<T, _nx, _nu, _np>::~Predictor() {
  }  

}


#endif
