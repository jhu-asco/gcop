#include <limits>
#include "uuvmanifold.h"
#include "se3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

UuvManifold::UuvManifold() : Manifold()
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
  //  SE3::Instance().log(dg, gi*xb.g);
  SE3::Instance().cayinv(dg, gi*xb.g);
  
  v.head<6>() = dg;
  v.tail<6>() = xb.v - xa.v;
}
 

void UuvManifold::Retract(UuvState &xb, 
                             const UuvState &xa,
                             const Vector12d &v) 
{
  Matrix4d dg;
  
  SE3::Instance().cay(dg, v.head<6>());
  
  xb.g = xa.g*dg;
  xb.v = xa.v + v.tail<6>();
}


void UuvManifold::dtau(Matrix12d &M, const Vector12d &v)
{
  M.setIdentity();
  Matrix6d D;
  SE3::Instance().dcay(D, v.head<6>());
  M.topLeftCorner<6,6>() = D;
}


void UuvManifold::Adtau(Matrix12d &M, const Vector12d &v)
{
  M.setIdentity();
  Matrix4d g;
  SE3::Instance().cay(g, v.head<6>());
  Matrix6d A;
  SE3::Instance().Ad(A, g);
  M.topLeftCorner<6,6>() = A;
}
