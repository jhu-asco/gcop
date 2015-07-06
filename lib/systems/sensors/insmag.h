#ifndef GCOP_INSMAG_H
#define GCOP_INSMAG_H

#include <Eigen/Dense>
#include "sensor.h"
#include "ins.h"
#include "rn.h"

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
  template<int _nu = 6, int _np = Dynamic>
    class InsMag : public Sensor<InsState, 15, _nu, _np, Matrix<double, 3, 1>, 3> {
    
  public:  
  
  typedef Matrix<double, _nu, 1> Vectorcd;

  typedef Matrix<double, _np, 1> Vectormd;


  typedef Matrix<double, 3, 1> Vectorzd;
  typedef Matrix<double, 3, 3> Matrixz3d;
  typedef Matrix<double, 3, _nu> Matrixzcd;
  typedef Matrix<double, 3, 15> Matrixz15d;
  typedef Matrix<double, 3, _np> Matrixzmd;
  
  InsMag();
  
  bool operator()(Vectorzd &y, double t, const InsState &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrixz15d *dydx = 0, Matrixzcd *dydu = 0,
                  Matrixzmd *dydp = 0) {
    

    y.template head<3>() = x.R.transpose()*m0;
    
    if (dydx){
      dydx->setZero();
      Matrix3d A;
      
      SO3::Instance().hat(A, y.template head<3>());
      dydx->template bottomLeftCorner<3,3>() = A;
    }
    
    return true;    
  }
  
  Vector3d m0;   ///< magnetometer reference
  
  double srm;   ///< magnetometer measurement stdev
  };
  
  
  template<int _nu, int _np>
    InsMag<_nu, _np>::InsMag() : 
    Sensor<InsState, 15, _nu, _np, Matrix<double, 3, 1>, 3>(Rn<3>::Instance()) {
    
    m0 << 1, 0, 0;
    
    srm = 0.01;       ///< magnetometer measurement stdev

    this->R.template bottomRightCorner<3,3>().diagonal().setConstant(srm*srm);
  }
}

#endif
