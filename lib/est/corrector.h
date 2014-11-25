#ifndef GCOP_CORRECTOR_H
#define GCOP_CORRECTOR_H

#include "sensor.h"
#include "manifold.h"

namespace gcop {
  
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic, 
    typename Tz = VectorXd,
    int _nz = Dynamic> class Corrector {
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
   * Basic Corrector
   * @param X sensor input manifold
   * @param sensor sensor 
   */
  Corrector(Manifold<T, _nx>  &X,
            Sensor<T, _nx, _nu, _np, Tz, _nz> &sensor);
  
  virtual ~Corrector();
    
  /**
   * Correction step. 
   * @param xb new belief state
   * @param t time
   * @param xa previous belief state
   * @param u control inputs
   * @param z measurement
   * @param p parameters (optional)
   * @param cov whether to update the covariance as well (true by default)
   * @return true if success
   */
  virtual bool Correct(T& xb, double t, const T &xa,
                       const Vectorcd &u, const Tz &z, 
                       const Vectormd *p = 0, bool cov = true) = 0;
    
  Manifold<T, _nx>  &X;             ///< system
  Sensor<T, _nx, _nu, _np, Tz, _nz> &sensor;  ///< sensor
    
  };
  
  
  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    Corrector<T, _nx, _nu, _np, Tz, _nz>::Corrector( Manifold<T, _nx>  &X, 
                                                     Sensor<T, _nx, _nu, _np, Tz, _nz> &sensor) : 
    X(X), sensor(sensor) {
  }  

  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    Corrector<T, _nx, _nu, _np, Tz, _nz>::~Corrector() {
  }  

}


#endif
