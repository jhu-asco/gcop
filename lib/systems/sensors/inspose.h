#ifndef GCOP_INSPOSE_H
#define GCOP_INSPOSE_H

#include <Eigen/Dense>
#include "sensor.h"
#include "posemanifold.h"

namespace gcop {
  
  using namespace Eigen;

  /**
   * Pose sensor: measures R and p directly 
   *
   * Subclasses should provide implementation for the 
   * sensor function ()
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template<int _nu = 6, int _np = Dynamic>
  class InsPose : public Sensor<InsState, 15, _nu, _np, PoseState, 6> {
    
  public:  

    typedef Matrix<double, _nu, 1> Vectorcd;
    typedef Matrix<double, _np, 1> Vectormd;

    typedef Matrix<double, 6, 15> Matrix6x15d;
    typedef Matrix<double, 6, _nu> Matrixrcd;
    typedef Matrix<double, 6, _np> Matrixrmd;
  
  InsPose();

  bool operator()(PoseState &y, double t, const InsState &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrix6x15d *dydx = 0, Matrixrcd *dydu = 0,
                  Matrixrmd *dydp = 0) {
    y.R = x.R;
    y.p = x.p;
    
    if (dydx){
      dydx->setZero();
      dydx->block<3,3>(0,0).setIdentity();
      dydx->block<3,3>(3,9).setIdentity();
    }
    
    return true;    
  }
  
  };
  
  template<int _nu, int _np> 
  InsPose<_nu, _np>::InsPose() : 
  Sensor<InsState, 15, _nu, _np, PoseState, 6>(PoseManifold::Instance()) {    
    
    this->R(0,0) = .001;
    this->R(1,1) = .001;
    this->R(2,2) = .001;
    this->R(3,3) = .001;
    this->R(4,4) = .001;
    this->R(5,5) = .001;
  }
}

#endif
