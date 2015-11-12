#include <limits>
#include "posemanifold.h"
#include "so3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

PoseManifold::PoseManifold() : Manifold()
{
}


PoseManifold& PoseManifold::Instance() 
{
  static PoseManifold instance;
  return instance;
}

void PoseManifold::Lift(Vector6d &v, 
                       const PoseState &xa,
                       const PoseState &xb) 
{
  Vector3d eR;
  //  SO3::Instance().log(eR, xa.R.transpose()*xb.R);
  SO3::Instance().cayinv(eR, xa.R.transpose()*xb.R);
  
  v.head<3>() = eR;
  v.segment<3>(3) = xb.p - xa.p;
}
 

void PoseManifold::Retract(PoseState &xb, 
                          const PoseState &xa,
                          const Vector6d &v) 
{
  Matrix3d dR;
  
  //  SO3::Instance().exp(dR, v.head<3>());
  SO3::Instance().cay(dR, v.head<3>());
  
  xb.R = xa.R*dR;
  xb.p = xa.p + v.segment<3>(3);
}


void PoseManifold::dtau(Matrix6d &M, const Vector6d &v)
{
  M.setIdentity();
  Matrix3d dR;
  SO3::Instance().dexp(dR, v.head<3>());
  M.topLeftCorner<3,3>() = dR;
}

void PoseManifold::Adtau(Matrix6d &M, const Vector6d &v)
{
  M.setIdentity();
  Matrix3d g;
  SO3::Instance().exp(g, v.head<3>());
  Matrix3d A;
  SO3::Instance().Ad(A, g);
  M.topLeftCorner<3,3>() = A;
}
