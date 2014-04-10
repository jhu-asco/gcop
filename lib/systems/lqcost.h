#ifndef GCOP_LQCOST_H
#define GCOP_LQCOST_H

#include <limits>
#include "cost.h"
#include <iostream>
#include "rn.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  template <typename T, int _n = Dynamic, int _c = Dynamic> class LqCost :
    public Cost<T, Matrix<double, _c, 1>, _n, _c>
    {
    public:
      
    typedef Matrix<double, _n, 1> Vectornd;
    typedef Matrix<double, _c, 1> Vectorcd;
    typedef Matrix<double, _n, _n> Matrixnd;
    typedef Matrix<double, _n, _c> Matrixncd;
    typedef Matrix<double, _c, _n> Matrixcnd;
    typedef Matrix<double, _c, _c> Matrixcd;
    
    /**
     * Linear-quadratic cost on a general manifold. Use this constructor for dynamic-size
     * control problem, i.e. LqCost<T>(X, U, tf, xf, ...)
     * @param X state space
     * @param U control space
     * @param tf final time
     * @param xf desired final state
     * @param diag whether the Q, R, and Qf matrices are diagonal?
     */
    LqCost(Manifold<T, _n> &X, Rn<_c> &U, double tf, const T &xf, bool diag = true);

    /**
     * Linear-quadratic cost on a general manifold. Use this constructor for fixed-size
     * control problem, i.e. LqCost<T, n, c>(X, tf, xf, ...)
     * @param tf final time
     * @param xf desired final state
     * @param diag whether the Q, R, and Qf matrices are diagonal?
     */
    LqCost(Manifold<T, _n> &X, double tf, const T &xf, bool diag = true);
    
    virtual double L(double t, const T& x, const Vectorcd& u, double h,
                     Vectornd *Lx = 0, Matrixnd* Lxx = 0,
                     Vectorcd *Lu = 0, Matrixcd* Luu = 0,
                     Matrixncd *Lxu = 0);

    /**
     * Set optional desired reference trajectory and controls
     * @param xds desired state trajectory (optional)
     * @param uds desired control sequence (optional)
     */
    void SetReference(const vector<T> *xds, const vector<Vectorcd> *uds);                          
    
    const T &xf; ///< reference to a desired final state
    
    bool diag;   ///< are the Q, Qf, and R matrices diagonal (true by default)
    
    Matrixnd Q;       ///< state matrix Q
    Matrixnd Qf;      ///< final state matrix Qf
    Matrixcd R;       ///< control matrix R

    const vector<T> *xds;         ///< optional reference trajectory 
    const vector<Vectorcd> *uds;  ///< optional reference control 

    protected:

    Vectornd dx;      ///< state error (as a tangent vector)
    Vectorcd du;      ///< control error

    };
  
  template <typename T, int _n, int _c>
    LqCost<T, _n, _c>::LqCost(Manifold<T, _n> &X, Rn<_c> &U, double tf, const T &xf, bool diag) : 
    Cost<T, Matrix<double, _c, 1>, _n, _c>(X, U, tf), xf(xf), diag(diag) {
    
    if (_n == Dynamic || _c == Dynamic) {
      Q.resize(X.n, X.n);
      Qf.resize(X.n, X.n);
      R.resize(U.n, U.n);
      dx.resize(X.n);
      du.resize(U.n);
    }

    //    Q.setZero();
    Q.setIdentity();
    Qf.setIdentity();
    R.setIdentity();

    xds = 0;
    uds = 0;
  }

  template <typename T, int _n, int _c>
    LqCost<T, _n, _c>::LqCost(Manifold<T, _n> &X, double tf, const T &xf, bool diag) : 
    Cost<T, Matrix<double, _c, 1>, _n, _c>(X, Rn<_c>::Instance(), tf), xf(xf), diag(diag) {
    
    assert(_c != Dynamic); 
    
    if (_n == Dynamic) {
      Q.resize(X.n, X.n);
      Qf.resize(X.n, X.n);
      dx.resize(X.n);
    }
    
    Q.setZero();
    Qf.setIdentity();
    R.setIdentity();
    
    xds = 0;
    uds = 0;
  }

  template <typename T, int _n, int _c>  
    void LqCost<T, _n, _c>::SetReference(const vector<T> *xds, const vector<Vectorcd> *uds) {
    this->xds = xds;
    this->uds = uds;
  }
  
  template <typename T, int _n, int _c> 
    double LqCost<T, _n, _c>::L(double t, const T &x, const Matrix<double, _c, 1> &u,
                                double h,
                                Matrix<double, _n, 1> *Lx, Matrix<double, _n, _n> *Lxx,
                                Matrix<double, _c, 1> *Lu, Matrix<double, _c, _c> *Luu,
                                Matrix<double, _n, _c> *Lxu) {

    int k = (int)(t/h);
    
    if (xds) {
      assert(k < xds->size());
      this->X.Lift(dx, (*xds)[k], x); // difference (on a vector space we have dx = x - xf)
    } else {
      this->X.Lift(dx, xf, x); // difference (on a vector space we have dx = x - xf)      
    }
    assert(!std::isnan(dx[0]));
    
    // check if final state
    if (t > this->tf - 1e-10) {
      if (Lx)
        if (diag)
          *Lx = Qf.diagonal().cwiseProduct(dx);
        else
          *Lx = Qf*dx;
      
      if (Lxx)
        *Lxx = Qf;
      
      if (Lu)
        Lu->setZero();
      if (Luu)
        Luu->setZero();
      if (Lxu)
        Lxu->setZero();    
      
      if (diag)
        return dx.dot(Qf.diagonal().cwiseProduct(dx))/2;
      else
        return dx.dot(Qf*dx)/2;
      
    } else {
      
      
      if (uds) {
        assert(k < uds->size());
        du = u - (*uds)[k];
      } else {
        du = u;
      }

      if (Lx)
        if (diag)
          *Lx = Q.diagonal().cwiseProduct(h*dx);
        else
          *Lx = Q*(h*dx);

      if (Lxx) 
        *Lxx = h*Q;
      
      if (Lu)
        if (diag)
          *Lu = R.diagonal().cwiseProduct(h*du);
        else
          *Lu = R*(h*du);        

      if (Luu)
        *Luu = h*R;

      if (Lxu)
        Lxu->setZero();
      
      if (diag)
        return h*(dx.dot(Q.diagonal().cwiseProduct(dx))/2 + du.dot(R.diagonal().cwiseProduct(du)))/2;
      else
        return h*(dx.dot(Q*dx) + du.dot(R*du))/2;
    }
    return 0;
  }  
}


#endif
