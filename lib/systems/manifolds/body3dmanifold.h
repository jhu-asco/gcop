// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_BODY3DMANIFOLD_H
#define GCOP_BODY3DMANIFOLD_H

#include "manifold.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 9, 1> Vector9d;
  typedef Matrix<double, 12, 1> Vector12d;
  typedef Matrix<double, 12, 12> Matrix12d;

  /**
   * Body3dState defined on SO(3)\times R^3 \times R^3 \times R^3
   *
   * Author: Marin Kobilarov
   */
  class Body3dState {
  public:
    Matrix3d R;    ///< rotation matrix
    Vector3d p;    ///< position
    Vector3d w;    ///< angular velocity
    Vector3d v;    ///< translational velocity

    void Clear() {
      R.setIdentity();
      p.setZero();
      w.setZero();
      v.setZero();
    }
  };
  
  class Body3dManifold : public Manifold<Body3dState, 12> {
    
  public:
    static Body3dManifold& Instance();

    void Lift(Vector12d &v,
              const Body3dState &xa,
              const Body3dState &xb);      

    void Retract(Body3dState &xb, 
                 const Body3dState &xa,
                 const Vector12d &v);

    void dtau(Matrix12d &M, const Vector12d &v);

    void dtauinv(Matrix12d &M, const Vector12d &v);

    void Adtau(Matrix12d &M, const Vector12d &v);

    bool useCay;   ///< whether to use the Cayley map instead of exponential

  private:
    Body3dManifold();
  };  
}


#endif
