#ifndef GCOP_KALMANPREDICTOR_H
#define GCOP_KALMANPREDICTOR_H

#include "predictor.h"

namespace gcop {
  
  template <typename T = VectorXd, 
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic> class KalmanPredictor : public Predictor<T, _nx, _nu, _np> {
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

  //  typedef pair<T, Matrixnd > Belief; ///< belief state 

  KalmanPredictor(System<T, _nx, _nu, _np> &sys);
  
  virtual ~KalmanPredictor();
    
  virtual bool Predict(T& xb, double t, const T &xa,
                       const Vectorcd &u, double h, 
                       const Vectormd *p = 0, bool cov = true);
  
  Matrixnd A; ///< process Jacobian
  Matrixnd Q; ///< discrete-time process noise covariance  
  };
  
  
  template <typename T, int _nx, int _nu, int _np> 
    KalmanPredictor<T, _nx, _nu, _np>::KalmanPredictor(System<T, _nx, _nu, _np>  &sys) : 
    Predictor<T, _nx, _nu, _np>(sys) {

    if (_nx == Dynamic) {
      A.resize(sys.X.n, sys.X.n);
      Q.resize(sys.X.n, sys.X.n);
    }
    A.setZero();
    Q.setZero();
  }
  
  template <typename T, int _nx, int _nu, int _np> 
    KalmanPredictor<T, _nx, _nu, _np>::~KalmanPredictor() {
  }  
  
  template <typename T, int _nx, int _nu, int _np> 
    bool KalmanPredictor<T, _nx, _nu, _np>::Predict(T& xb, double t, const T &xa,
                                                    const Vectorcd &u, double h, 
                                                    const Vectormd *p, bool cov) {
    this->sys.Step(xb, t, xa, u, h, p, &A);
    if (cov) {
      this->sys.Noise(Q, t, xa, u, h, p);
      xb.P = A*xa.P*A.transpose() + Q;
    }
    return true;
  }
}


#endif
