#ifndef GCOP_QROTOR_H
#define GCOP_QROTOR_H

#include "body3d.h"

namespace gcop {
  
   /**
   * A quadrotor dynamical model.
   *
   * The state is
   * \f$ \bf x = (R, x, \omega, v) \f$ where \f$ (R,x)\in SO(3)\times\mathbb{R}^3\f$ is the pose, 
   * and the controls are \f$ \bm u = (u_1,u_2,u_3,u_4)\f$ correspond to torques 
   * around the body and a vertical lift force. 
   *
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  class Qrotor : public Body3d<4> {
  public:
    Qrotor();    
    
    double l;  ///< distance from center of mass to each rotor
    double r;  ///< propeller radius
    
    double kt; ///< thrust gain
    double km; ///< moment gain
    
    double mm; ///< motor coefficient
  };
}

#endif
