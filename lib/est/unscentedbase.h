#ifndef GCOP_UNSCENTEDBASE_H
#define GCOP_UNSCENTEDBASE_H

#include "manifold.h"

namespace gcop {

  using namespace std;
  
  template <typename T = VectorXd, 
    int _nx = Dynamic> class UnscentedBase {
  public:
  
  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nx, _nx> Matrixnd;

  UnscentedBase(Manifold<T, _nx> &X);
  
  virtual ~UnscentedBase();

  /**
   * Compute 2*L+1 sigma points given mean x and spread A
   * @param Xs a vector of 2*L+1 states
   * @param x mean state
   * @param A the square root of the covariance 
   */
  void Points(vector<T> &xs,
              const T& x);
  
  Manifold<T, _nx> &X; ///< manifold

  double a;      ///< UKF \f$\alpha\f$-param (a=1e-3 by default)
  double k;      ///< UKF \f$\kappa\f$-param (k=0 by default)
  double b;      ///< UKF \f$\beta\f$-param (b=2 by default)
  
  protected:
  int L;         ///< number of sigma points (this is set internally to 2*model.nr + 1
  double l;      ///< the variable \f$ l = a^2(L+k)-L \f$
  
  vector<T> Xs;   ///< state sigma points
  
  vector<double> Ws;  ///< points weights
  vector<double> Wc;  ///< cov update weights

  Matrixnd A;   ///< cholesky factor
  };
  

  
  template <typename T, int _nx> 
    UnscentedBase<T, _nx>::UnscentedBase(Manifold<T, _nx> &X) :
    X(X),
    Xs(2*X.n + 1),
    Ws(2*X.n + 1),
    Wc(2*X.n + 1) {
    
    this->L = X.n;
    this->a = .01;
    // this->a = 1;
    this->k = 0;
    this->b = 2;     
    //    this->b = 0;
    this->l = a*a*(L+k)-L;
    
    for (int i = 0; i < 2*L + 1; ++i) {      
      if (i == 0) {
        Ws[0] = l/(L+l);
        Wc[0] = l/(L+l) + (1-a*a+b);
      } else {
        Ws[i] = 1/(2*(L+l));
        Wc[i] = 1/(2*(L+l));
      }
    }
          
    if (_nx == Dynamic) {
      A.resize(X.n, X.n);
    }
  }
  
  template <typename T, int _nx> 
    UnscentedBase<T, _nx>::~UnscentedBase() {
  }  
  template <typename T, int _nx> 
    void UnscentedBase<T, _nx>::Points(vector<T> &Xs,
                                       const T &x) {
    Xs[0] = x;
    
    A = x.P.llt().matrixL().transpose();
    bool pd = true;
    
    if (!pd) {
      cout << "[W] UKF::Predict: Cholesky failed!" << endl;
    }
    
    //    cout << "A=" << A << endl;
    
    for (int i = 0; i < L; ++i) {
      Vectornd dx = sqrt(L+l)*A.col(i);
      Vectornd dxm = -sqrt(L+l)*A.col(i);
      
      this->X.Retract(Xs[i + 1], x, dx);
      /*
      cout << "#i=" << i << endl;
      cout << "x.R=" << x.R << endl;
      cout << "dx=" << dx.transpose() << endl;
      cout << "R=" << Xs[i+1].R << endl;
      */
      this->X.Retract(Xs[i + 1 + L], x, dxm);
    }
  }
}


#endif
