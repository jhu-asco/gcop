// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_BODY2D_H
#define GCOP_BODY2D_H

#include "system.h"
#include <utility>
#include "body2dforce.h"
#include "body2dmanifold.h"
#include "se2.h"
#include "rn.h"


namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 6, 3> Matrix63d;
  typedef Matrix<double, 3, 6> Matrix36d;
  
  typedef pair<Matrix3d, Vector3d> Body2dState;

   /**
   * A "geometric" planar rigid body system. The control system
   * is implemented as an evolution on SE(2) as opposed to a vector space.
   *
   * The state is
   * \f$ \bm x = (g, v) \f$ where \f$ g\in SE(2)\f$ is the pose, 
   * and \f$ v\in \mathbb{R}^3\f$ is the velocity
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  template<int c = 3>
  class Body2d : public System<Body2dState, 6, c> {
  public:

  typedef Matrix<double, c, 1> Vectorcd;
  typedef Matrix<double, 6, c> Matrix6cd;
  typedef Matrix<double, 3, c> Matrix3cd;

  Body2d(Body2dForce<c> *f = 0);
  
  double Step(Body2dState &xb, double t, const Body2dState &xa,
              const Vectorcd &u, double h, const VectorXd *p = 0,
              Matrix6d *A = 0, Matrix6cd *B = 0, Matrix<double, 6, Dynamic> *C = 0);
  
  Vector2d d; ///< body dimensions
  Vector3d I; ///< inertia components
  
  Body2dForce<c> *force;  ///< force from controls and external disturbances
  
  };
  
  template<int c>
    Body2d<c>::Body2d(Body2dForce<c> *force) : 
    System<Body2dState, 6, c>(Body2dManifold::Instance()), 
    d(1, .3), I(.3,10,10), force(force)
  {
  }
  
  template<int c>
     double Body2d<c>::Step(Body2dState& xb, double t, const Body2dState& xa,
                            const Vectorcd& u, double h, const VectorXd* p,
                            Matrix6d *A, Matrix6cd *B, Matrix<double, 6, Dynamic> *C) {
    Matrix36d fx;  // force Jacobian with respect to x
    Matrix3cd fu;    // force Jacobian with respect to u
    Matrix<double, 3, Dynamic> fp;  // force Jacobian with respect to p
    
    if (p && C)
    fp.resize(3, p->size());

    if (force) {
      Vector3d f;
      force->Set(f, xa, t, u, h, p,
                 A ? &fx : 0,
                 B ? &fu : 0,
                 C ? &fp : 0);
      xb.second = xa.second + h*f;
    } else {
      xb.second = xa.second;
    }
    
    Vector3d hxib = h*xb.second;
    
    SE2 &se2 = SE2::Instance();
    
    Matrix3d dg;
    
    se2.cay(dg, hxib);
    xb.first = xa.first*dg;
    
    if (A) {
      se2.cay(dg, -hxib);
      Matrix3d Adm;
      se2.Ad(Adm, dg);
      A->topLeftCorner<3,3>() = Adm;
      se2.dcay(dg, -hxib);    
      A->topRightCorner<3,3>() = h*dg;        
      if (force) {
        A->bottomRows<3>() = h*fx;
      }
      A->bottomRightCorner<3,3>() += Matrix3d::Identity();
    }
    
    if (B) {
      // A should've been called too
      assert(A);
      
      B->topRows(3) = (h*h)*(dg*fu);
      //  B->bottomRows(3) = (Vector3d(h, h, h).cwiseQuotient(I)).asDiagonal()*fu;
      B->bottomRows(3) = h*fu;
      
    }
    
    if (C) {
      if (p) {
        assert(C->cols() == p->size() && C->rows() == 6); 
      }
      if (force) {
        assert(A);
        C->topRows<3>() = (h*h)*(dg*fp);
        C->bottomRows<3>() = (Vector3d(h, h, h).cwiseQuotient(I)).asDiagonal()*fp;
      } else {
        C->setZero();
      }
    }
    
  return 1;
}


}


#endif
