// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <limits>
#include "pose3dmanifold.h"
#include "so3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Pose3dManifold::Pose3dManifold() : Manifold(), useCay(false)
{
}


Pose3dManifold& Pose3dManifold::Instance() 
{
  static Pose3dManifold instance;
  return instance;
}

void Pose3dManifold::Lift(Vector6d &v, 
                          const Pose3d &xa,
                          const Pose3d &xb) 
{
  const Matrix3d &Ra = xa.R;  
  const Matrix3d &Rb = xb.R;
  
  Vector3d eR;

  if (useCay)
    SO3::Instance().cayinv(eR, Ra.transpose()*Rb);
  else
    SO3::Instance().log(eR, Ra.transpose()*Rb);
  
  v.head<3>() = eR;
  v.tail<3>() = xb.p - xa.p;
}
 

void Pose3dManifold::Retract(Pose3d &xb, 
                             const Pose3d &xa,
                             const Vector6d &v) 
{
  Matrix3d dR;
  
  if (useCay)
    SO3::Instance().cay(dR, v.head<3>());
  else
    SO3::Instance().exp(dR, v.head<3>());
  
  xb.R = xa.R*dR;
  xb.p = xa.p + v.tail<3>();
}


void Pose3dManifold::dtau(Matrix6d &M, const Vector6d &v)
{
  M.setIdentity();
  Matrix3d dR;
  
  if (useCay)
    SO3::Instance().dcay(dR, v.head<3>());
  else
    SO3::Instance().dexp(dR, v.head<3>());
  
  M.topLeftCorner<3,3>() = dR;
}

void Pose3dManifold::dtauinv(Matrix6d &M, const Vector6d &v)
{
  M.setIdentity();
  Matrix3d dR;
  if (useCay)
    SO3::Instance().dcayinv(dR, v.head<3>());
  else
    SO3::Instance().dexpinv(dR, v.head<3>());
  
  M.topLeftCorner<3,3>() = dR;
}


void Pose3dManifold::Adtau(Matrix6d &M, const Vector6d &v)
{
  M.setIdentity();
  Matrix3d g;

  if (useCay)
    SO3::Instance().cay(g, v.head<3>());
  else
    SO3::Instance().exp(g, v.head<3>());
    
  Matrix3d A;
  SO3::Instance().Ad(A, g);
  M.topLeftCorner<3,3>() = A;
}
