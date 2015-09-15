#ifndef GCOP_INSPOSE3D_H
#define GCOP_INSPOSE3D_H

#include <Eigen/Dense>
#include "sensor.h"
#include "ins.h"
#include "pose3dmanifold.h"

namespace gcop {
  
  using namespace Eigen;
  
  /**
   * pose3d sensor for the Ins system
   *
   * Subclasses should provide implementation for the 
   * sensor function ()
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template<int _nu, int _np = Dynamic>
    class InsPose3d : public Sensor<InsState, 15, _nu, _np,
    Pose3d, 6> {
    
  public:  
    
    typedef Matrix<double, _nu, 1> Vectorcd;    
    typedef Matrix<double, _np, 1> Vectormd;
        
    typedef Matrix<double, 6, 15> Matrix6x15d;
    typedef Matrix<double, 6, 3> Matrix6x3d;
    typedef Matrix<double, 6, _nu> Matrix6cd;
    typedef Matrix<double, 6, _np> Matrix6md;
    
  InsPose3d();
  
  bool operator()(Pose3d &y, double t, const InsState &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrix6x15d *dydx = 0, Matrix6cd *dydu = 0,
                  Matrix6md *dydp = 0) {    
    y.R = x.R;
    y.p = x.p;
    
    if (dydx){
      dydx->block<3,3>(0,0).setIdentity();
      dydx->block<3,3>(9,3).setIdentity();
    }
    
    return true;    
  }
  
  Vector3d sp;  ///< position measurement standard deviation
  Vector3d sR;  ///< rotation measurement standard deviation

  
  };
  
  
  template<int _nu, int _np>
    InsPose3d<_nu, _np>::InsPose3d() :
    Sensor<InsState, 15, _nu, _np, Pose3d, 6>(Pose3dManifold::Instance()) {
    
    sp.setConstant(.1);
    sR.setConstant(.1);
    
    this->R.template topLeftCorner<3,3>().diagonal() = sp.cwiseProduct(sp);
    this->R.template bottomRightCorner<3,3>().diagonal() = sR.cwiseProduct(sp);
  }
}

#endif //GCOP_INSPOSE3D_H
