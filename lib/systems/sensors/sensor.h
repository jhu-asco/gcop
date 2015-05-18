#ifndef GCOP_SENSOR_H
#define GCOP_SENSOR_H

#include <Eigen/Dense>
#include <type_traits>
#include <iostream>
#include "manifold.h"

namespace gcop {
  
  using namespace Eigen;
    
  /**
   * General sensor model 
   *
   * Subclasses should provide implementation for the 
   * sensor function ()
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <typename T = VectorXd, 
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic, 
    typename Tz = T,
    int _nz = _nx> class Sensor {
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
    
  typedef Matrix<double, _nz, 1> Vectorrd;
  typedef Matrix<double, _nz, _nz> Matrixrd;
  typedef Matrix<double, _nz, _nx> Matrixrnd;
  typedef Matrix<double, _nz, _nu> Matrixrcd;
  typedef Matrix<double, _nz, _np> Matrixrmd;

  Sensor(Manifold<Tz, _nz> &Z);

  /**
   * Measurement functions
   * @param y deterministic output (i.e. the actual random measurement is z = y + v, where the noise v is Normal(0,var) )
   * @param t current time
   * @param x current state
   * @param u control inputs (often not used)
   * @param p static parameters  (optional)
   * @param dydx jacobian w.r.t. x   (optional)
   * @param dydu jacobian w.r.t. u   (optional)
   * @param dydp jacobian w.r.t. p   (optional)
   * @return true if all variables are feasible
   */
  virtual bool operator()(Tz &y, double t, const T &x, const Vectorcd &u,
                          const Vectormd *p = 0, 
                          Matrixrnd *dydx = 0, Matrixrcd *dydu = 0,
                          Matrixrmd *dydp = 0);
  
  Manifold<Tz, _nz> &Z;        ///< measurement manifold  
  Matrixrd R;                  ///< sensor noise covariance
  
  };


  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    Sensor<T, _nx, _nu, _np, Tz, _nz>::Sensor(Manifold<Tz, _nz> &Z) : Z(Z) {
    if (_nz == Dynamic) {
      R.resize(Z.n, Z.n);
    }
    R.setIdentity();
  }
  
  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    bool Sensor<T, _nx, _nu, _np, Tz, _nz>::operator()(Tz &y, 
                                                       double t, const T &x, const Vectorcd &u,
                                                       const Vectormd *p, 
                                                       Matrixrnd *dydx, Matrixrcd *dydu,
                                                       Matrixrmd *dydp) {
    if(std::is_same<Tz, T>::value)//Default Implementation if T and Tz are same
    {
        y = (Tz)(*((Tz *)(&x)));//Hacky way of copying pointer over
        //y = (Tz)x;
        return true;
    }
    std::cout << "[W] Sensor::(): unimplemented! Subclasses should override." << std::endl;
    return false;
  }
}

#endif

