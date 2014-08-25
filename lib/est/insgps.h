#ifndef GCOP_INSGPS_H
#define GCOP_INSGPS_H

#include <Eigen/Dense>
#include "sensor.h"
#include "ins.h"

namespace gcop {
  
  using namespace Eigen;

  //  typedef Insgps<6> FullInsgps;
  //  typedef Insgps<3> AccInsgps;
  
  /**
   * General sensor model 
   *
   * Subclasses should provide implementation for the 
   * sensor function ()
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class InsGps : public Sensor<InsState, 15, 6, Dynamic, Vector3d, 3> {
    
  public:  
  
  typedef Matrix<double, 3, 6> Matrix36d;
  typedef Matrix<double, 3, 15> Matrix3x15d;
  typedef Matrix<double, 3, Dynamic> Matrix3Xd;
  
  InsGps();

  bool operator()(Vector3d &y, double t, const InsState &x, const Vector6d &u,
                  const VectorXd *p = 0, 
                  Matrix3x15d *dydx = 0, Matrix36d *dydu = 0,
                  Matrix3Xd *dydp = 0) {
    
    //    y.template head<3>() = x.first.transpose()*a0 - ba;
    y = x.p;
    
    if (dydx){
      dydx->block<3,3>(0,9).setIdentity();
    }
    
    return true;    
  }
  
  double sxy;   ///< x-y position standard deviation
  double sz;    ///< altitute stdev
  };
  
  
  InsGps::InsGps() : 
  Sensor<InsState, 15, 6, Dynamic, Vector3d, 3>(Rn<3>::Instance()) {    
    
    sxy = 3;
    sz = 10;
    
    this->R(0,0) = sxy*sxy;
    this->R(1,1) = sxy*sxy;
    this->R(2,2) = sz*sz;
  }
}

#endif
