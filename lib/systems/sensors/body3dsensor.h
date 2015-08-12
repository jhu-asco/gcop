#ifndef GCOP_BODY3DSENSOR_H
#define GCOP_BODY3DSENSOR_H

#include <Eigen/Dense>
#include "sensor.h"
#include "body3d.h"
#include "pose3dmanifold.h"

namespace gcop {
  
  using namespace Eigen;

  //  typedef Body3dsensor<6> FullBody3dsensor;
  //  typedef Body3dsensor<3> AccBody3dsensor;
  
  /**
   * Body3d sensor model 
   *
   * Subclasses should provide implementation for the 
   * sensor function ()
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template<int _nu, int _np = Dynamic>
    class Body3dSensor : public Sensor<Body3dState, 12, _nu, _np, 
    Pose3d, 6> {
    
  public:  
    
    typedef Matrix<double, _nu, 1> Vectorcd;    
    typedef Matrix<double, _np, 1> Vectormd;
        
    typedef Matrix<double, 6, 12> Matrix6x12d;
    typedef Matrix<double, 6, 3> Matrix6x3d;
    typedef Matrix<double, 6, _nu> Matrix6cd;
    typedef Matrix<double, 6, _np> Matrix6md;
    
  Body3dSensor();
  
  bool operator()(Pose3d &y, double t, const Body3dState &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrix6x12d *dydx = 0, Matrix6cd *dydu = 0,
                  Matrix6md *dydp = 0) {    
    y.R = x.R;
    y.p = x.p;
    
    if (dydx){
      dydx->setIdentity();
    }
    
    return true;    
  }
  
  Vector3d sp;  ///< position measurement standard deviation
  Vector3d sR;  ///< rotation measurement standard deviation

  
  };
  
  
  template<int _nu, int _np>
    Body3dSensor<_nu, _np>::Body3dSensor() : 
    Sensor<Body3dState, 12, _nu, _np, Pose3d, 6>(Pose3dManifold::Instance()) {
    
    sp.setConstant(.1);
    sR.setConstant(.1);
    
    this->R.template topLeftCorner<3,3>().diagonal() = sp.cwiseProduct(sp);
    this->R.template bottomRightCorner<3,3>().diagonal() = sR.cwiseProduct(sp);
  }
}

#endif
