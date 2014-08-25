#ifndef GCOP_UNSCENTEDCORRECTOR_H
#define GCOP_UNSCENTEDCORRECTOR_H

#include "corrector.h"
#include "unscentedbase.h"

namespace gcop {
  
  using namespace std;
  
  template <typename T = VectorXd, 
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic, 
    typename Tz = VectorXd,
    int _nz = Dynamic> class UnscentedCorrector : public Corrector<T, _nx, _nu, _np, Tz, _nz>, public UnscentedBase<T, _nx> {
  public:
  
  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;

  typedef Matrix<double, _nx, _nx> Matrixnd;
  typedef Matrix<double, _nx, _nu> Matrixncd;
  typedef Matrix<double, _nx, _nz> Matrixnrd;
  typedef Matrix<double, _nu, _nx> Matrixcnd;
  typedef Matrix<double, _nu, _nu> Matrixcd;
  
  typedef Matrix<double, _np, _np> Matrixmd;
  typedef Matrix<double, _nx, _np> Matrixnmd;
  typedef Matrix<double, _np, _nx> Matrixmnd;

  typedef Matrix<double, _nz, 1> Vectorrd;
  typedef Matrix<double, _nz, _nz> Matrixrd;
  typedef Matrix<double, _nz, _nx> Matrixrnd;
  typedef Matrix<double, _nz, _nu> Matrixrcd;
  typedef Matrix<double, _nz, _np> Matrixrmd;

  UnscentedCorrector(System<T, _nx, _nu, _np> &sys,
                     Sensor<T, _nx, _nu, _np, Tz, _nz> &sensor);
  
  virtual ~UnscentedCorrector();
  
  virtual bool Correct(T& xb, double t, const T &xa, 
                       const Vectorcd &u, const Tz &z, 
                       const Vectormd *p = 0, bool cov = true);    

  protected:

  vector<Tz> Zs;  ///< measurement sigma points  
  Matrixrd Pzz;   ///< internally used covariances
  Matrixnrd Pxz;  ///< internally used covariances
  
  };
  
  
  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    UnscentedCorrector<T, _nx, _nu, _np, Tz, _nz>::UnscentedCorrector(System<T, _nx, _nu, _np>  &sys, 
                                                                      Sensor<T, _nx, _nu, _np, Tz, _nz> &sensor) : 
    Corrector<T, _nx, _nu, _np, Tz, _nz>(sys, sensor),
    UnscentedBase<T, _nx>(sys.X),
    Zs(2*sys.X.n + 1) {
    
    if (_nz == Dynamic) {     
      Pzz.resize(sensor.Z.n, sensor.Z.n);
    }      
    if (_nx == Dynamic || _nz == Dynamic) {     
      Pxz.resize(sys.X.n, sensor.Z.n);        
    }
  }
  
  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    UnscentedCorrector<T, _nx, _nu, _np, Tz, _nz>::~UnscentedCorrector() {
  }  
  

  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    bool UnscentedCorrector<T, _nx, _nu, _np, Tz, _nz>::Correct(T& xb, double t, const T &xa, 
                                                                const Vectorcd &u, const Tz &z, 
                                                                const Vectormd *p, bool cov) {
    
    this->Points(this->Xs, xa);

    //    VectorXd xm = VectorXd::Zero(x.size());
    Vectorrd zm;
    zm.setZero();

    Vectornd dx = Vectornd::Zero();
    vector<Vectornd> dxs(2*this->L+1);

    for (int i = 0; i < 2*this->L+1; ++i) {
      sensor(Zs[i], t, this->Xs[i], u, p);
      zm = zm + this->Ws[i]*Zs[i];
      //      dx = dx + Ws[i]*dxs[i];
      //cout << "Zs["<<i<<"]=" << Zs[i].transpose() << endl;
      //      cout << "Ws["<<i<<"]=" << Ws[i] << endl;
    }

    cout << "z(0)=" << Zs[0].transpose() <<endl;
    cout << "zm=" << zm.transpose() << endl;

    Pzz.setZero();
    Pxz.setZero();    

    //    T xm;
    //    this->sys.X.Retract(xm, xa, dx);

    for (int i = 0; i < 2*this->L+1; ++i) {
      Vectorrd dz = Zs[i] - zm;
      this->sys.X.Lift(dxs[i], xa, this->Xs[i]);      // variations from mean
      /*
      Vector3d rpy;
      SO3::Instance().g2q(rpy, Xs[i].R);
      cout << "rpy=" << rpy.transpose() << endl;

      cout << "dxs["<<i<<"]=" << dxs[i].transpose() << endl;
      */

      //      this->sys.X.Lift(dxs[i], xm, Xs[i]);
      Pzz = Pzz + this->Wc[i]*dz*dz.transpose();
      Pxz = Pxz + this->Wc[i]*dxs[i]*dz.transpose();
    }
    cout <<"Pzz=" << Pzz << endl;
    cout <<"Pxz=" << Pxz << endl;

    // Pxz is like P*H'

    //    cout << "z=" << z.transpose() << endl;

    Pzz = Pzz + this->sensor.R;

    Matrixnrd K = Pxz*Pzz.inverse();
    cout << "K=" << K << endl;
    
    dx = K*(z - zm);
    this->sys.X.Retract(xb, xa, dx);
    
    xb.P = xa.P - K*Pzz*K.transpose();
    //    xb.P = xa.P - K*Pxz.transpose();
  
    return true;
  }
}


#endif
