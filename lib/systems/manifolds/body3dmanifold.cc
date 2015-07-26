// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <limits>
#include "body3dmanifold.h"
#include "so3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Body3dManifold::Body3dManifold() : Manifold(), useCay(false)
{
}


Body3dManifold& Body3dManifold::Instance() 
{
  static Body3dManifold instance;
  return instance;
}

void Body3dManifold::Lift(Vector12d &v, 
                          const Body3dState &xa,
                          const Body3dState &xb) 
{
  const Matrix3d &Ra = xa.R;  
  const Matrix3d &Rb = xb.R;
  
  Vector3d eR;

  if (useCay)
    SO3::Instance().cayinv(eR, Ra.transpose()*Rb);
  else
    SO3::Instance().log(eR, Ra.transpose()*Rb);
  
  v.head<3>() = eR;
  v.segment<3>(3) = xb.p - xa.p;
  v.segment<3>(6) = xb.w - xa.w;
  v.tail<3>() = xb.v - xa.v;
}
 

void Body3dManifold::Retract(Body3dState &xb, 
                             const Body3dState &xa,
                             const Vector12d &v) 
{
  Matrix3d dR;
  
  if (useCay)
    SO3::Instance().cay(dR, v.head<3>());
  else
    SO3::Instance().exp(dR, v.head<3>());
  
  xb.R = xa.R*dR;
  xb.p = xa.p + v.segment<3>(3);
  xb.w = xa.w + v.segment<3>(6);
  xb.v = xa.v + v.tail<3>();
}


void Body3dManifold::dtau(Matrix12d &M, const Vector12d &v)
{
  M.setIdentity();
  Matrix3d dR;

  if (useCay)
    SO3::Instance().dcay(dR, v.head<3>());
  else
    SO3::Instance().dexp(dR, v.head<3>());

  M.topLeftCorner<3,3>() = dR;
}

void Body3dManifold::dtauinv(Matrix12d &M, const Vector12d &v)
{
  M.setIdentity();
  Matrix3d dR;
  if (useCay)
    SO3::Instance().dcayinv(dR, v.head<3>());
  else
    SO3::Instance().dexpinv(dR, v.head<3>());

  M.topLeftCorner<3,3>() = dR;
}


void Body3dManifold::Adtau(Matrix12d &M, const Vector12d &v)
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
