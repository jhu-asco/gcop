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
   * A hexrotor with two manipulators
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  class Airbot : public Mbs {
  public:
    Airbot();    

    void Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
               MatrixXd *A = 0, MatrixXd *B = 0);
    
    Hrotor hrotor;
    
  };
}

#endif
