#ifndef GCOP_UNSCENTEDPREDICTOR_H
#define GCOP_UNSCENTEDPREDICTOR_H

#include "predictor.h"
#include "unscentedbase.h"

namespace gcop {

  using namespace std;
  
  template <typename T = VectorXd, 
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic> class UnscentedPredictor : public Predictor<T, _nx, _nu, _np>, public UnscentedBase<T, _nx> {
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

  UnscentedPredictor(System<T, _nx, _nu, _np> &sys);
  
  virtual ~UnscentedPredictor();
  
  virtual bool Predict(T& xb, double t, const T &xa,
                       const Vectorcd &u, double h, 
                       const Vectormd *p = 0, bool cov = true);
  
  
  protected:

  vector<T> Xps;  ///< predicted state sigma points
  Matrixnd Q;     ///< discrete-time process noise covariance  
  
  };
  
  
  template <typename T, int _nx, int _nu, int _np> 
    UnscentedPredictor<T, _nx, _nu, _np>::UnscentedPredictor(System<T, _nx, _nu, _np>  &sys) : 
    Predictor<T, _nx, _nu, _np>(sys),
    UnscentedBase<T, _nx>(sys.X),
    Xps(2*sys.X.n + 1) {
    
    if (_nx == Dynamic) {
      Q.resize(sys.X.n, sys.X.n);
    }
    Q.setZero();
  }
  
  template <typename T, int _nx, int _nu, int _np> 
    UnscentedPredictor<T, _nx, _nu, _np>::~UnscentedPredictor() {
  }  
  
  template <typename T, int _nx, int _nu, int _np> 
    bool UnscentedPredictor<T, _nx, _nu, _np>::Predict(T& xb, double t, const T &xa,
                                                       const Vectorcd &u, double h, 
                                                       const Vectormd *p, bool cov) {
    
    this->Points(this->Xs, xa);
    
    // average difference
    Vectornd dx = Vectornd::Zero();
    vector<Vectornd> dxs(2*this->L+1);
    
    for (int i = 0; i < 2*this->L + 1; ++i) {
      this->sys.Step(Xps[i], t, this->Xs[i], u, h, p);

      /*
      Vector3d rpy;
      SO3::Instance().g2q(rpy, Xs[i].R);
      cout << "rpy=" << rpy.transpose() << endl;

      SO3::Instance().g2q(rpy, Xps[i].R);
      cout << "rpyp=" << rpy.transpose() << endl;
      */


      this->sys.X.Lift(dxs[i], Xps[0], Xps[i]);   
      dx = dx + this->Ws[i]*dxs[i];   ///< average in exponential coordinates
    }

    this->sys.X.Retract(xb, Xps[0], dx);
    xb.P.setZero(); // zero covariance
    
    for (int i = 0; i < 2*this->L+1; ++i) {
      this->sys.X.Lift(dxs[i], xb, Xps[i]);   ///< variation from propagated mean
      xb.P = xb.P + this->Wc[i]*dxs[i]*dxs[i].transpose();
    }

    this->sys.Noise(Q, t, xa, u, h, p);
    xb.P = xb.P + Q;

    return true;
  }  
}


#endif
