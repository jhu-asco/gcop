// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_BODY3DCOST_H
#define GCOP_BODY3DCOST_H

#include "body3d.h"
#include "lqcost.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  /**
   * Linear-quadratic cost for rigid body systems
   *
   * Author: Marin Kobilarov
   */
  template <int c = 6> class Body3dCost : public LqCost<Body3dState, 12, c, Dynamic, 6> {
    
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, c, c> Matrixcd;
    typedef Matrix<double, 6, c> Matrix6xcd;
    typedef Matrix<double, 12, c> Matrix12xcd;
    
  public:
    
    Body3dCost(Body3d<c> &sys, double tf, const Body3dState &xf, bool diag = true);
  };  
  
  template <int c> Body3dCost<c>::Body3dCost(Body3d<c> &sys, double tf, const Body3dState &xf, bool diag) : 
    LqCost<Body3dState, 12, c, Dynamic, 6>(sys, tf, xf, diag) {

    this->Qf(0,0) = .5;
    this->Qf(1,1) = .5;
    this->Qf(2,2) = .5;
    this->Qf(3,3) = 5;
    this->Qf(4,4) = 5;
    this->Qf(5,5) = 5;
    
    this->Qf(6,6) = .1;
    this->Qf(7,7) = .1;
    this->Qf(8,8) = .1;
    this->Qf(9,9) = 1;
    this->Qf(10,10) = 1;
    this->Qf(11,11) = 1;
    
    this->R.diagonal() = Matrix<double, c, 1>::Constant(.1);
  }
}


#endif
