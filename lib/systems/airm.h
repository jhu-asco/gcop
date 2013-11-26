#ifndef GCOP_AIRM_H
#define GCOP_AIRM_H

#include "mbs.h"
#include "hrotor.h"

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
  class Airm : public Mbs {
  public:
    Airm();    

    void Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
               MatrixXd *A = 0, MatrixXd *B = 0);
    
    Hrotor hrotor;
    
  };
}

#endif
