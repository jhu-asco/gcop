#ifndef GCOP_LQCOST_H
#define GCOP_LQCOST_H

#include <limits>
#include "lscost.h"
#include <iostream>
#include <type_traits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  template <typename T, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic,
    int _ng = Dynamic> class LqCost : public LsCost<T, _nx, _nu, _np, _ng, T> {
    public:

  typedef Matrix<double, _ng, 1> Vectorgd;
  typedef Matrix<double, _ng, _nx> Matrixgxd;
  typedef Matrix<double, _ng, _nu> Matrixgud;
  typedef Matrix<double, _ng, _np> Matrixgpd;

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;

  typedef Matrix<double, _nx, _nx> Matrixnd;
  typedef Matrix<double, _nx, _nu> Matrixncd;
  typedef Matrix<double, _nu, _nx> Matrixcnd;
  typedef Matrix<double, _nu, _nu> Matrixcd;
  
  //  typedef Matrix<double, Dynamic, 1> Vectormd;
  typedef Matrix<double, _np, _np> Matrixmd;
  typedef Matrix<double, _nx, _np> Matrixnmd;
  typedef Matrix<double, _np, _nx> Matrixmnd;
      
    /**
     * Linear-quadratic cost on a general manifold. Use this constructor for dynamic-size
     * control problem, i.e. LqCost<T>(X, U, tf, xf, ...)
     * @param X state space
     * @param U control space
     * @param tf final time
     * @param xf desired final state
     * @param diag whether the Q, R, and Qf matrices are diagonal?
     */
    LqCost(System<T, _nx, _nu, _np> &sys, double tf, const T &xf, bool diag = true);

    /**
     * Updates the square root gains, i.e. sets Qsqrt = Q.sqrt(), etc...
     *
     */
    void UpdateGains(); 
    
    virtual double L(double t, const T& x, const Vectorcd& u, double h,
                     const Vectormd *p = 0,
                     Vectornd *Lx = 0, Matrixnd* Lxx = 0,
                     Vectorcd *Lu = 0, Matrixcd* Luu = 0,
                     Matrixncd *Lxu = 0,
                     Vectormd *Lp = 0, Matrixmd *Lpp = 0,
                     Matrixmnd *Lpx = 0);

    bool Res(Vectorgd &g, 
             double t, const T &x, const Vectorcd &u, double h,
             const Vectormd *p = 0,
             Matrixgxd *dgdx = 0, Matrixgud *dgdu = 0,
             Matrixgpd *dgdp = 0);                 

    /**
     * Set optional desired reference trajectory and controls
     * @param xds desired state trajectory (optional)
     * @param uds desired control sequence (optional)
     */
    void SetReference(const vector<T> *xds, const vector<Vectorcd> *uds);        
    
    virtual bool SetContext(const T &c);

    const T *xf; ///< reference to a desired final state
    
    bool diag;   ///< are the Q, Qf, and R matrices diagonal (true by default)
    
    Matrixnd Q;       ///< state matrix Q
    Matrixnd Qf;      ///< final state matrix Qf
    Matrixcd R;       ///< control matrix R

    Matrixnd Qsqrt;       ///< state matrix Q sqrt
    Matrixnd Qfsqrt;      ///< final state matrix Qf sqrt
    Matrixcd Rsqrt;       ///< control matrix R sqrt

    const vector<T> *xds;         ///< optional reference trajectory 
    const vector<Vectorcd> *uds;  ///< optional reference control 

    protected:

    Vectornd dx;      ///< state error (as a tangent vector)
    Vectorcd du;      ///< control error
    Matrixnd D;       ///< used for computing jacobians

    };
  
  template <typename T, int _nx, int _nu, int _np, int _ng>
    LqCost<T, _nx, _nu, _np, _ng>::LqCost(System<T, _nx, _nu, _np> &sys, 
                                          double tf, const T &xf, bool diag) : 
    LsCost<T, _nx, _nu, _np, _ng, T>(sys, tf, sys.X.n + sys.U.n), xf(&xf), diag(diag) {
    
    if (_nx == Dynamic || _nu == Dynamic) {
      Q.resize(sys.X.n, sys.X.n);
      Qf.resize(sys.X.n, sys.X.n);
      R.resize(sys.U.n, sys.U.n);
      Qsqrt.resize(sys.X.n, sys.X.n);
      Qfsqrt.resize(sys.X.n, sys.X.n);
      Rsqrt.resize(sys.U.n, sys.U.n);
      dx.resize(sys.X.n);
      du.resize(sys.U.n);
      D.resize(sys.X.n, sys.X.n);
    }

    Q.setZero();
    //    Q.setIdentity();
    Qf.setIdentity();
    R.setIdentity();

    UpdateGains();

    D.setIdentity();

    xds = 0;
    uds = 0;
  }

  template <typename T, int _nx, int _nu, int _np, int _ng>
    bool LqCost<T, _nx, _nu, _np, _ng>::SetContext(const T& c) {
    this->xf = &c;
  }
  
  template <typename T, int _nx, int _nu, int _np, int _ng>
    void LqCost<T, _nx, _nu, _np, _ng>::UpdateGains() {
    Qsqrt = Q.array().sqrt();
    Qfsqrt = Qf.array().sqrt();
    Rsqrt = R.array().sqrt();
  }


  template <typename T, int _nx, int _nu, int _np, int _ng>  
    void LqCost<T, _nx, _nu, _np, _ng>::SetReference(const vector<T> *xds, const vector<Vectorcd> *uds) {
    this->xds = xds;
    this->uds = uds;
  }

  template <typename T, int _nx, int _nu, int _np, int _ng> 
    bool LqCost<T, _nx, _nu, _np, _ng>::Res(Vectorgd &g, 
                                            double t, const T &x, const Vectorcd &u, double h,
                                            const Vectormd *p,
                                            Matrixgxd *dgdx, Matrixgud *dgdu,
                                            Matrixgpd *dgdp) {
    assert(g.size() == this->sys.U.n + this->sys.X.n);
    
    int &nx = this->sys.X.n;
    int &nu = this->sys.U.n;
    
    // check if final state
    if (t > this->tf - 1e-10) {
      if (xds) {
        this->sys.X.Lift(dx, xds->back(), x); // difference (on a vector space we have dx = x - xf)
      } else {
        this->sys.X.Lift(dx, *xf, x); // difference (on a vector space we have dx = x - xf)     
      }
      if (std::isnan(dx[0])) {
        cout << "[W] LqCost::Res: dx is nan" << endl;
        return false;
      }

      if (diag)
        g.head(nx) = Qfsqrt.diagonal().cwiseProduct(dx);
      else
        g.head(nx) = Qfsqrt*dx;      

      g.tail(nu).setZero();

    } else {
      assert(h > 0);
      int k = round(t/h);
      if (xds) {
        assert(k < xds->size());
        this->sys.X.Lift(dx, (*xds)[k], x); // difference (on a vector space we have dx = x - xf)
      } else {
        this->sys.X.Lift(dx, *xf, x); // difference (on a vector space we have dx = x - xf)     
      }
      if (std::isnan(dx[0])) {
        cout << "[W] LqCost::Res: dx is nan" << endl;
        return false;
      }

      if (diag)
        g.head(nx) = Qsqrt.diagonal().cwiseProduct(sqrt(h)*dx);
      else
        g.head(nx) = Qsqrt*(sqrt(h)*dx);

      if (uds) {
        assert(k < uds->size());
        du = sqrt(h)*(u - (*uds)[k]);
      } else {
        du = sqrt(h)*u;
      }
      
      if (diag)
        g.tail(nu) = Rsqrt.diagonal().cwiseProduct(du);
      else
        g.tail(nu) = Rsqrt*(du);
    }
    return true;
  }
  
  
  template <typename T, int _nx, int _nu, int _np, int _ng> 
    double LqCost<T, _nx, _nu, _np, _ng>::L(double t, const T &x, const Matrix<double, _nu, 1> &u,
                                            double h,
                                            const Matrix<double, _np, 1> *p,
                                            Matrix<double, _nx, 1> *Lx, Matrix<double, _nx, _nx> *Lxx,
                                            Matrix<double, _nu, 1> *Lu, Matrix<double, _nu, _nu> *Luu,
                                            Matrix<double, _nx, _nu> *Lxu,
                                            Matrix<double, _np, 1> *Lp, Matrix<double, _np, _np> *Lpp,
                                            Matrix<double, _np, _nx> *Lpx) {
    
    int k = (h > 0.0) ? round(t/h) : 0 ;
    
        // check if final state
    if (t > this->tf - 1e-10) {

      if (xds) {
        this->sys.X.Lift(dx, xds->back(), x); // difference (on a vector space we have dx = x - xf)
      } else {
        this->sys.X.Lift(dx, *xf, x); // difference (on a vector space we have dx = x - xf)      
      }
      if (std::isnan(dx[0])) {
        cout << "[W] LqCost::L: dx is nan" << endl;
        return std::numeric_limits<double>::max();
      }    

      if (Lx) {
        if (diag)
          *Lx = Qf.diagonal().cwiseProduct(dx);
        else
          *Lx = Qf*dx;
        
        if (!std::is_same<T, Matrix<double, _nx, 1> >::value) {
          this->sys.X.dtauinv(D, -dx);
          *Lx = D.transpose()*(*Lx);
        }

      }
      if (Lxx) {
        if (!std::is_same<T, Matrix<double, _nx, 1> >::value) {
          assert(Lx); // Lx should also be requested, so D is already computed
          // sys.X.dtauinv(D, -dx);
          *Lxx = D.transpose()*Qf*D;
        } else {
          *Lxx = Qf;      
        }
      }

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

      if (xds) {
        assert(k < xds->size());
        this->sys.X.Lift(dx, (*xds)[k], x); // difference (on a vector space we have dx = x - xf)
      } else {
        this->sys.X.Lift(dx, *xf, x); // difference (on a vector space we have dx = x - xf)      
      }
      if (std::isnan(dx[0])) {
        cout << "[W] LqCost::L: dx is nan" << endl;
        return std::numeric_limits<double>::max();
      }     
      
      if (uds) {
        assert(k < uds->size());
        du = u - (*uds)[k];
      } else {
        du = u;
      }

      if (Lx) {
        if (diag)
          *Lx = Q.diagonal().cwiseProduct(h*dx);
        else
          *Lx = Q*(h*dx);

        // if not Euclidean space
        if (!std::is_same<T, Matrix<double, _nx, 1> >::value) {
          this->sys.X.dtauinv(D, -dx);
          *Lx = D.transpose()*(*Lx);
        }

      }
      if (Lxx) {
        if (!std::is_same<T, Matrix<double, _nx, 1> >::value) {
          assert(Lx); // Lx should also be requested, so D is already computed
          //          sys.X.dtauinv(D, -dx);
          *Lxx = h*(D.transpose()*Q*D);
        } else {          
          *Lxx = h*Q;
        }
      }

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
        return h*(dx.dot(Q.diagonal().cwiseProduct(dx)) + du.dot(R.diagonal().cwiseProduct(du)))/2;
      else
        return h*(dx.dot(Q*dx) + du.dot(R*du))/2;
    }
    return 0;
  }  
}


#endif
