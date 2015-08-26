// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef GCOP_UUVMANIFOLD_H
#define GCOP_UUVMANIFOLD_H

#include "manifold.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 9, 1> Vector9d;
  typedef Matrix<double, 12, 1> Vector12d;
  typedef Matrix<double, 12, 12> Matrix12d;
  

  /**
   *  State on SE(3)\times R^6
   *
   */
  struct UuvState {
    UuvState() {g.setIdentity(); v.setZero(); }
    virtual ~UuvState() {};

    Matrix4d g;  ///< SE(3) configuration
    Vector6d v;  ///< body-fixed velocities
  };
  

  /**
   * Left-trivialized TSE(3) manifold
   *
   * Author: Marin Kobilarov
   */
  class UuvManifold : public Manifold<UuvState, 12> {
    
  public:
    static UuvManifold& Instance();

    void Lift(Vector12d &v,
              const UuvState &xa,
              const UuvState &xb);      

    void Retract(UuvState &xb, 
                 const UuvState &xa,
                 const Vector12d &v);

    void dtau(Matrix12d &M, const Vector12d &v);

    void dtauinv(Matrix12d &M, const Vector12d &v);

    void Adtau(Matrix12d &M, const Vector12d &v);

    char tauType;                  ///< type of tau map (default is TAU_EXP)

    static const int TAU_EXP = 1;  ///< exponential map
    static const int TAU_CAY = 2;  ///< Cayley map

  private:
    UuvManifold();
  };  
}


#endif
