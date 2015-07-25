#include <limits>
#include "uuvmanifold.h"
#include "se3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

UuvManifold::UuvManifold() : Manifold(), tauType(TAU_EXP)
{
}


UuvManifold& UuvManifold::Instance() 
{
  static UuvManifold instance;
  return instance;
}

void UuvManifold::Lift(Vector12d &v, 
                       const UuvState &xa,
                       const UuvState &xb) 
{
  Vector6d dg;
  Matrix4d gi;
  SE3::Instance().inv(gi, xa.g);
  
  switch(tauType) {
  case TAU_EXP:
    SE3::Instance().log(dg, gi*xb.g);
    break;
  case TAU_CAY:
    SE3::Instance().cayinv(dg, gi*xb.g);
    break;
  default:
    break;
  }
  v.head<6>() = dg;
  v.tail<6>() = xb.v - xa.v;
}
 

void UuvManifold::Retract(UuvState &xb, 
                             const UuvState &xa,
                             const Vector12d &v) 
{
  Matrix4d dg;
  
  switch(tauType) {
  case TAU_EXP:
    SE3::Instance().exp(dg, v.head<6>());
    break;
  case TAU_CAY:
    SE3::Instance().cay(dg, v.head<6>());
    break;
  default:
    break;
  }
  
  xb.g = xa.g*dg;
  xb.v = xa.v + v.tail<6>();
}


void UuvManifold::dtau(Matrix12d &M, const Vector12d &v)
{
  M.setIdentity();
  Matrix6d D;
  switch(tauType) {
  case TAU_EXP:
    SE3::Instance().dexp(D, v.head<6>());
    break;
  case TAU_CAY:
    SE3::Instance().dcay(D, v.head<6>());
    break;
  default:
    break;
  }
  M.topLeftCorner<6,6>() = D;
}

void UuvManifold::dtauinv(Matrix12d &M, const Vector12d &v)
{
  M.setIdentity();
  Matrix6d D;
  switch(tauType) {
  case TAU_EXP:
    SE3::Instance().dexpinv(D, v.head<6>());
    break;
  case TAU_CAY:
    SE3::Instance().dcayinv(D, v.head<6>());
    break;
  default:
    break;
  }
  M.topLeftCorner<6,6>() = D;
}


void UuvManifold::Adtau(Matrix12d &M, const Vector12d &v)
{
  M.setIdentity();
  Matrix4d g;
  switch(tauType) {
  case TAU_EXP:
    SE3::Instance().exp(g, v.head<6>());
    break;
  case TAU_CAY:
    SE3::Instance().cay(g, v.head<6>());
    break;
  default:
    break;
  }
  Matrix6d A;
  SE3::Instance().Ad(A, g);
  M.topLeftCorner<6,6>() = A;
}
