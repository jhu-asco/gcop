#ifndef GCOP_INSIMU_H
#define GCOP_INSIMU_H

#include <Eigen/Dense>
#include "sensor.h"
#include "ins.h"
#include "rn.h"

namespace gcop {
  
  using namespace Eigen;

  //  typedef InsImu<6> FullInsImu;
  //  typedef InsImu<3> AccInsImu;
  
  /**
   * General sensor model 
   *
   * Subclasses should provide implementation for the 
   * sensor function ()
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template<int _nz = 6, int _nu = 6, int _np = Dynamic>
    class InsImu : public Sensor<InsState, 15, _nu, _np, Matrix<double, _nz, 1>, _nz> {
    
  public:  
  
  typedef Matrix<double, _nu, 1> Vectorcd;

  typedef Matrix<double, _np, 1> Vectormd;


  typedef Matrix<double, _nz, 1> Vectorzd;
  typedef Matrix<double, _nz, 3> Matrixz3d;
  typedef Matrix<double, _nz, _nu> Matrixzcd;
  typedef Matrix<double, _nz, 15> Matrixz15d;
  typedef Matrix<double, _nz, _np> Matrixzmd;
  
  InsImu();
  
  bool operator()(Vectorzd &y, double t, const InsState &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrixz15d *dydx = 0, Matrixzcd *dydu = 0,
                  Matrixzmd *dydp = 0) {
    
    //    y.template head<3>() = x.first.transpose()*a0 - ba;
    y.template head<3>() = x.R.transpose()*a0 - x.ba;

    if (_nz == 6)
      y.template tail<3>() = x.R.transpose()*m0;
    
    if (dydx){
      dydx->setZero();
      Matrix3d A;
      SO3::Instance().hat(A, y.template head<3>());
      dydx->template topLeftCorner<3,3>() = A;
      dydx->template block<3,3>(0,6).diagonal().setConstant(-1);
      
      // for magnetometer measurements
      if (_nz == 6) {
        SO3::Instance().hat(A, y.template tail<3>());
        dydx->template bottomLeftCorner<3,3>() = A;
      }
    }
    
    return true;    
  }
  
  Vector3d a0;   ///< acceleration reference
  Vector3d m0;   ///< magnetometer reference
  
  double sra;   ///< acceleration measurement stdev
  double srm;   ///< magnetometer measurement stdev
  };
  
  
  template<int _nz, int _nu, int _np>
    InsImu<_nz, _nu, _np>::InsImu() : 
    Sensor<InsState, 15, _nu, _np, Matrix<double, _nz, 1>, _nz>(Rn<_nz>::Instance()) {
    
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
