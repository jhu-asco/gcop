// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_UUV_H
#define GCOP_UUV_H

#include "system.h"
#include "uuvmanifold.h"
#include "se3.h"
#include <limits>
#include <iostream>
#include <utility>
#include "function.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 12, 12> Matrix12d;
  
  /**
   * An underwater/sufrace vehicle system modeled as a single rigid body on SE(3)
   *
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <int c = 6> class Uuv : public System<UuvState, 12, c> {
    
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, c, c> Matrixcd;
    typedef Matrix<double, 6, c> Matrix6xcd;
    typedef Matrix<double, 12, c> Matrix12xcd;
    typedef Matrix<double, 12, Dynamic> Matrix12Xd;
    
  public:
    
    Uuv();

    virtual ~Uuv();
    
    double Step(UuvState &xb, double t, const UuvState &xa, 
                const Vectorcd &u, double h, const VectorXd *p = 0,
                Matrix12d *A = 0, Matrix12xcd *B = 0, Matrix12Xd *C = 0);

    double SympStep(double t, UuvState &xb, const UuvState &xa, 
                    const Vectorcd &u, double h, const VectorXd *p = 0,
                    Matrix12d *A = 0, Matrix12xcd *B = 0, Matrix12Xd *C = 0);
    
    double EulerStep(double t, UuvState &xb, const UuvState &xa, 
                     const Vectorcd &u, double h, const VectorXd *p = 0,
                     Matrix12d *A = 0, Matrix12xcd *B = 0, Matrix12Xd *C = 0);    

    static void Compute(Matrix6d &I, double m, const Vector3d &ds);
    
    Matrix6d I;    ///< 6x6 inertia matrix

    Vector6d d;    ///< 6x1 decoupled quadratic drag terms
    
    Vector3d b;    ///< center of boyancy

    Vector3d fp;   ///< constant position force in spatial frame (e.g. due to gravity)
    
    Matrix6xcd B; ///< control input transformation

    bool symp;     ///< symplectic?
		
    string name;   ///< Unique name of the body

    Vector3d ds;   ///< body dimensions (currently used only for visualization)

  };  
  

  
  template <int c> 
    Uuv<c>::Uuv() : 
    System<UuvState, 12, c>(UuvManifold::Instance()),
    I(Matrix6d::Identity()),
    d(Vector6d::Zero()),
    b(0,0,0),
    fp(0,0,0),
    B(Matrix<double, 6, c>::Identity()),
    symp(true),
    ds(1.5,.8,.5) {    
    Compute(I, 50, ds);
  }
  
  template<int c>  Uuv<c>::~Uuv()
    {
    }  

  template <int c>
    double Uuv<c>::Step(UuvState &xb, double t, const UuvState &xa, 
                           const Matrix<double, c, 1> &u, double h, const VectorXd *p,
                           Matrix12d *A, Matrix<double, 12, c> *B, Matrix12Xd *C) {
    if (symp)
      return SympStep(t, xb, xa, u, h, p, A, B, C);
    else
      return EulerStep(t, xb, xa, u, h, p, A, B, C);      
  }
  
  template <int c>
    double Uuv<c>::SympStep(double t, UuvState &xb, const UuvState &xa, 
                            const Matrix<double, c, 1> &u, double h, const VectorXd *p,
                            Matrix12d *A, Matrix<double, 12, c> *B, Matrix12Xd *C) {
    
    SE3 &se3 = SE3::Instance();
    
    Matrix6d Da, Db;
    
    // initialize wb
    xb.v = xa.v;

    // at this point the following can be iterated more than once

    se3.tln(Da, -h*xa.v);
    se3.tln(Db, h*xb.v);

    Vector6d fd = d.array()*(xa.v.array().abs()*xa.v.array() + xb.v.array().abs()*xb.v.array())/2;  // quadratic drag

    const Matrix3d &Ra = xa.g.topLeftCorner<3,3>();

    Vector6d fe;
    fe.head<3>() = b.cross(Ra.row(2));  // boyancy force
    fe.tail<3>() = Ra.transpose()*fp;  // gravity force

    // dynamic residual
    Vector6d e = Db.transpose()*I*xb.v - Da.transpose()*I*xa.v - h*(fd + fe + this->B*u);

    Matrix6d adv;
    Matrix6d admu;
    se3.ad(adv, xb.v); 
    se3.adt(admu, I*xb.v); 

    Vector6d dfd = 2*d.array()*xb.v.array().abs();
    Matrix6d De = I - h/2*(adv.transpose()*I + admu + Matrix6d(dfd.asDiagonal()));
    
    xb.v = xb.v - De.inverse()*e;    

    Matrix4d dg;
    se3.cay(dg, xb.v);
    xb.g = xa.g*dg;
    
    return 1;
  }


  /////////////////////////////////////////////////

  template <int c>
    double Uuv<c>::EulerStep(double t, UuvState &xb, const UuvState &xa, 
                                const Matrix<double, c, 1> &u, double h, const VectorXd *p,
                                Matrix12d *A, Matrix<double, 12, c> *B, Matrix12Xd *C) {
    return 1;
  }

  template <int c> 
    void Uuv<c>::Compute(Matrix6d &I, double m, const Vector3d &ds) {
    I.setZero();
    I(0,0) = m*(ds[1]*ds[1] + ds[2]*ds[2])/3;
    I(1,1) = m*(ds[0]*ds[0] + ds[2]*ds[2])/3;
    I(2,2) = m*(ds[0]*ds[0] + ds[1]*ds[1])/3;
    I(3,3) = m; 
    I(4,4) = m;
    I(5,5) = m;
  }


}
#endif
