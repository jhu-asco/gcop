#ifndef GCOP_FILTER_H
#define GCOP_FILTER_H

#include "system.h"
#include "sensor.h"

namespace gcop {
  
  template <typename T = VectorXd, 
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic, 
    typename Tz = VectorXd,
    int _nz = Dynamic> class Filter {
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
   * Basic Filter
   * @param model system model
   */
  Filter(System<T, _nx, _nu, _np>  &sys,
         Sensor<T, _nx, _nu, _np, Tz, _nz> &sensor);
  
  virtual ~Filter();
    
  /**
   * Prediction step. 
   * @param u control 
   * @param cov whether to update the covariance as well (true by default)
   * @return predicted state
   */
  virtual bool Predict(T& xb, double t, const T &xa, 
                       const Vectorcd &u, double h, 
                       const Vectormd *p = 0, bool cov = true) = 0;
  
  /**
   * Correction step. The base implementation does not change the state.
   * @param z measurement
   * @param cov whether to update the covariance as well (true by default)
   * @return corrected state
   */
  virtual bool Update(T& xb, double t, const T &xa,
                      const Vectorcd &u, const Tz &z, 
                      const Vectormd *p = 0, bool cov = true) = 0;
  
  /**
   * Set measurement Chi Gate
   * @param chi squared chi bound
   */
  void SetChiGate(double chiGate) { this->chiGate = chiGate; };
  
  /**
   * Get the preset measurement Chi Gate value
   * @return chi squared chi bound
   */
  double GetChiGate() const {return chiGate; }
  
  /**
   * Returns the internally computed chi-value during the call to Correct
   * @return the chi-value of the most recent measurement
   */
  double GetChi() const {return chi; }
  
  System<T, _nx, _nu, _np>  &sys;             ///< system
  Sensor<T, _nx, _nu, _np, Tz, _nz> &sensor;  ///< sensor
    
  double chi;             ///< computed chi-estimate (computed only if chiGate > 0)
  double chiGate;         ///< chi-gate (default: zero or negative to ignore)
  
  protected:

  };
  
  
  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    Filter<T, _nx, _nu, _np, Tz, _nz>::Filter( System<T, _nx, _nu, _np>  &sys, 
                                               Sensor<T, _nx, _nu, _np, Tz, _nz> &sensor) : 
    sys(sys), sensor(sensor) {
  }  

  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    Filter<T, _nx, _nu, _np, Tz, _nz>::~Filter() {
  }  

}


#endif
