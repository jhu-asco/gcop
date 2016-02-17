// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_QROTOR_ID_MODEL_COST_H
#define GCOP_QROTOR_ID_MODEL_COST_H

#include "qrotoridmodel.h"
#include "lqcost.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  /**
   * Linear-quadratic cost for Quadrotor system
   *
   * Author: Marin Kobilarov
   */
  class QRotorIDModelCost : public LqCost<QRotorIDState, 15, 4, 13, 19> {
    
    typedef Matrix<double, 4, 1> Vectorcd;
    typedef Matrix<double, 4, 4> Matrixcd;
    typedef Matrix<double, 6, 4> Matrix6xcd;
    typedef Matrix<double, 12, 4> Matrix12xcd;
    
  public:
    
    QRotorIDModelCost(QRotorIDModel &sys, double tf, const QRotorIDState &xf, bool diag = true)
        : LqCost<QRotorIDState, 15, 4, 13, 19>(sys, tf, xf, diag) {

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
    this->Qf(12,12) = 0;
    this->Qf(13,13) = 0;
    this->Qf(14,14) = 0;

    this->R.diagonal() = Matrix<double, 4, 1>::Constant(.1);
    this->UpdateGains();
  }
  };
}


#endif
