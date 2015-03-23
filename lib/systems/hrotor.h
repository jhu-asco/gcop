#ifndef GCOP_HROTOR_H
#define GCOP_HROTOR_H

#include "body3d.h"

namespace gcop {
  
   /**
   * A hexrotor dynamical model.
   *
   * The state is
   * \f$ \bf x = (R, x, \omega, v) \f$ where \f$ (R,x)\in SO(3)\times\mathbb{R}^3\f$ is the pose, 
   * and the controls are \f$ \bm u = (u_1,u_2,u_3,u_4)\f$ correspond to torques 
   * around the body and a vertical lift force.
   *
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  class Hrotor : public Body3d<4> {
  public:
    Hrotor();
    void StateAndControlsToFlat(VectorXd &y, const Body3dState &x,
               const Vector4d &u);    
    void FlatToStateAndControls(Body3dState &x, Vector4d &u,
               const std::vector<VectorXd> &y);
    
    double l;  ///< distance from center of mass to each rotor
    double r;  ///< propeller radius
    
    double kt; ///< thrust gain
    double km; ///< moment gain
    
    double mm; ///< motor coefficient
  };
}

#endif
