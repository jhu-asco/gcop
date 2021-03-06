#ifndef GCOP_IMUSENSOR_H
#define GCOP_IMUSENSOR_H

#include <Eigen/Dense>
#include "sensor.h"
#include "imu.h"

namespace gcop {
  
  using namespace Eigen;
  
  /**
   * An IMU Sensor model. Initialize with nz=3 for no magnetometer, otherwise
   * set nz=6 (default) for acceleration + magnetometer measurements.
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template<int _nz = 6, int _nu = 3, int _np = Dynamic>
    class ImuSensor : public Sensor<ImuState, 9, _nu, _np, Matrix<double, _nz, 1>, _nz> {
    
  public:  
  
  typedef Matrix<double, _nz, 1> Vectorzd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;

  typedef Matrix<double, _nz, _nu> Matrixzcd;
  typedef Matrix<double, _nz, 6> Matrixz6d;
  typedef Matrix<double, _nz, 9> Matrixz9d;
  typedef Matrix<double, _nz, _np> Matrixzmd;
  
  ImuSensor();

  bool operator()(Vectorzd &y, double t, const ImuState &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrixz9d *dydx = 0, Matrixzcd *dydu = 0,
                  Matrixzmd *dydp = 0) {
    
    y.template head<3>() = x.R.transpose()*a0 - x.ba;

    if (_nz == 6)
      y.template tail<3>() = x.R.transpose()*m0;
    
    if (dydx){
      dydx->setZero();
      Matrix3d A;
      SO3::Instance().hat(A, y.template head<3>());
      dydx->template topLeftCorner<3,3>() = A;
      dydx->template topRightCorner<3,3>().diagonal().setConstant(-1);
      
      // for magnetometer measurements
      if (_nz == 6) {
        SO3::Instance().hat(A, y.template tail<3>());
        dydx->template bottomLeftCorner<3,3>() = A;
      }
    }
    
    return true;    
  }
  
  Vector3d a0;   ///< acceleration reference (e.g. Earth gravity)
  Vector3d m0;   ///< magnetometer reference (e.g. magnetic North)
  
  double sra;   ///< acceleration measurement stdev
  double srm;   ///< magnetometer measurement stdev
  };
  
  
  template<int _nz, int _nu, int _np>
    ImuSensor<_nz, _nu, _np>::ImuSensor() : 
    Sensor<ImuState, 9, _nu, _np, Matrix<double, _nz, 1>, _nz>(Rn<_nz>::Instance()) {
    
    a0 << 0, 0, -9.81; 
    m0 << 1, 0, 0;     
    
    sra = 0.01;       ///< acceleration measurement stdev
    srm = 0.01;       ///< magnetometer measurement stdev

    this->R.template topLeftCorner<3,3>().diagonal().setConstant(sra*sra);
    if (_nz == 6)
      this->R.template bottomRightCorner<3,3>().diagonal().setConstant(srm*srm);

  }
}

#endif

