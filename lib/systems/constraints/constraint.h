// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_CONSTRAINT_H
#define GCOP_CONSTRAINT_H

#include <Eigen/Dense>
#include <assert.h>

namespace gcop {
  
  using namespace Eigen;  
  
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic,
    int _ng = Dynamic> class Constraint {
    
  public:
  
  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;
  typedef Matrix<double, _ng, 1> Vectorgd;

  typedef Matrix<double, _ng, _nx> Matrixgnd;
  typedef Matrix<double, _ng, _nu> Matrixgcd;
  typedef Matrix<double, _ng, _np> Matrixgmd;

  typedef Matrix<double, _nx, _nx> Matrixnd;
  typedef Matrix<double, _nx, _nu> Matrixncd;
  typedef Matrix<double, _nu, _nx> Matrixcnd;
  typedef Matrix<double, _nu, _nu> Matrixcd;
  
  //  typedef Matrix<double, Dynamic, 1> Vectormd;
  typedef Matrix<double, _np, _np> Matrixmd;
  typedef Matrix<double, _nx, _np> Matrixnmd;
  typedef Matrix<double, _np, _nx> Matrixmnd;
    
  Constraint(int ng = 0) : ng(_ng != Dynamic ? _ng : ng) {
    assert(this->ng > 0);    
  }

  // constraint g(t,x,u,p)<=0 
  virtual bool operator()(Vectorgd &g,
                          double t, const T &x, const Vectorcd &u,
                          const Vectormd *p = 0, 
                          Matrixgnd *dgdx = 0, Matrixgcd *dgdu = 0,
                          Matrixgmd *dgdp = 0) = 0;
  
  // constraint g(t,x)<=0 
  virtual bool operator()(Vectorgd &g,
                          double t, const T &x, 
                          Matrixgnd *dgdx = 0) {
    return (*this)(g, t, x, ub, 0, dgdx);
  }

  int ng;        ///< constraint dimension
  
  private:
  Vectorcd ub;   ///< blank control 
  };
}

#endif

