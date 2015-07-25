#ifndef GCOP_INSGPS_H
#define GCOP_INSGPS_H

#include <Eigen/Dense>
#include "sensor.h"
#include "ins.h"

namespace gcop {
  
  using namespace Eigen;

  /**
   * GPS sensor 
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template<int _nu = 6, int _np = Dynamic>
  class InsGps : public Sensor<InsState, 15, _nu, _np, Vector3d, 3> {
    
  public:  

    typedef Matrix<double, _nu, 1> Vectorcd;
    typedef Matrix<double, _np, 1> Vectormd;

    typedef Matrix<double, 3, 15> Matrix3x15d;
    typedef Matrix<double, 3, _nu> Matrixrcd;
    typedef Matrix<double, 3, _np> Matrixrmd;
  
  InsGps();

  bool operator()(Vector3d &y, double t, const InsState &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrix3x15d *dydx = 0, Matrixrcd *dydu = 0,
                  Matrixrmd *dydp = 0) {
    
    y = x.p;
    
    if (dydx){
      dydx->block<3,3>(0,9).setIdentity();
    }
    
    return true;    
  }
  
  double sxy;   ///< x-y position standard deviation
  double sz;    ///< altitute stdev
  };
  
  template<int _nu, int _np> 
  InsGps<_nu, _np>::InsGps() : 
  Sensor<InsState, 15, _nu, _np, Vector3d, 3>(Rn<3>::Instance()) {    
    
    sxy = 3;
    sz = 10;
    
    this->R(0,0) = sxy*sxy;
    this->R(1,1) = sxy*sxy;
    this->R(2,2) = sz*sz;
  }
}

#endif
