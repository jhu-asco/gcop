// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_SYSTEM_H
#define GCOP_SYSTEM_H

#include <Eigen/Dense>
#include <vector>

namespace gcop {
  
  using namespace Eigen;

  /**
   * A container for generic trajectory (this could be discrete or continuous) but in any case
   * we contain a discrete-time representation for computational/visualization purposes
   *
   * Author: Marin Kobilarov (c) marin(at)jhu.edu
   */
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic> class Trajectory {
  public:
  typedef Matrix<double, _nu, 1> Vectorud;
  typedef Matrix<double, _np, 1> Vectorpd;
  
  Trajectory(const System<T, _nx, _nu, _np> &sys, int N);
  
  const System<T, _nx, _nu, _np> &sys;  ///< system 
  vector<double> ts;                    ///< sequence of (N+1) times
  vector<T> xs;                         ///< sequence of (N+1) states
  vector<Vectorud> us;                  ///< sequence of (N) control inputs (regarded as paremetrizing the control signal along the i-th segment, for i=0,...,N-1) e.g. constant control during the interval \f$[t_k,t_{k+1}]\f$
  Vectorpd *p;                          ///< system parameters (optional, set to zero by default)
  
  };

  template <typename T, int _nx, int _nu, int _np> 
    Trajectory<T, _nx, _nu, _np>::Trajectory(const System<T, _nx, _nu, _np> &sys, int N) : 
    sys(sys), ts(N+1), xs(N+1), us(N), p(0)
    {      
    }  
}

#endif
