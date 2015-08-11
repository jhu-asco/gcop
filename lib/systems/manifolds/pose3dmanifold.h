// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_POSE3DMANIFOLD_H
#define GCOP_POSE3DMANIFOLD_H

#include "manifold.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;

  /**
   * Pose3d defined on SO(3)\times R^3
   *
   * Author: Marin Kobilarov
   */
  class Pose3d {
  public:
    Matrix3d R;    ///< rotation matrix
    Vector3d p;    ///< position

    Matrix6d P;   ///< covariance

    void Clear() {
      R.setIdentity();
      p.setZero();
      P.setIdentity();
    }

    void resize(int n) {};
  };
  
  class Pose3dManifold : public Manifold<Pose3d, 6> {
    
  public:
    static Pose3dManifold& Instance();

    void Lift(Vector6d &v,
              const Pose3d &xa,
              const Pose3d &xb);      
    
    void Retract(Pose3d &xb, 
                 const Pose3d &xa,
                 const Vector6d &v);
    
    void dtau(Matrix6d &M, const Vector6d &v);
    
    void dtauinv(Matrix6d &M, const Vector6d &v);
    
    void Adtau(Matrix6d &M, const Vector6d &v);
    
    bool useCay;   ///< whether to use the Cayley map instead of exponential
    
  private:
    Pose3dManifold();
  };  
}


#endif
