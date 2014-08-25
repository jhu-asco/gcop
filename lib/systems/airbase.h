// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_AIRBOT_H
#define GCOP_AIRBOT_H

#include "mbs.h"
#include "hrotor.h"


namespace gcop {
  
   /**
   * A quadrotor based mbs model.
   *
   * The state is
   * \f$ \bf x = (R, x, \omega, v) \f$ where \f$ (R,x)\in SO(3)\times\mathbb{R}^3\f$ is the pose, 
   * and the controls are \f$ \bm u = (u_1,u_2,u_3,u_4)\f$ correspond to torques 
   * around the body and a vertical lift force. 
   *
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  class Airbase : public Mbs {
  public:
    Airbase(int nb, int j); ///< nb number of bodies ///<j number of joints 

    void Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
               MatrixXd *A = 0, MatrixXd *B = 0);
    
    
  };
}

#endif
