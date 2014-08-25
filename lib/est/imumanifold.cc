#include <limits>
#include "imumanifold.h"
#include "so3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

ImuManifold::ImuManifold() : Manifold()
{
}


ImuManifold& ImuManifold::Instance() 
{
  static ImuManifold instance;
  return instance;
}

void ImuManifold::Lift(Vector9d &v, 
                       const ImuState &xa,
                       const ImuState &xb) 
{
  Vector3d eR;
  //  SO3::Instance().log(eR, xa.R.transpose()*xb.R);
  SO3::Instance().cayinv(eR, xa.R.transpose()*xb.R);
  
  v.head<3>() = eR;
  v.segment<3>(3) = xb.bg - xa.bg;
  v.tail<3>() = xb.ba - xa.ba;
}
 

void ImuManifold::Retract(ImuState &xb, 
                          const ImuState &xa,
                          const Vector9d &v) 
{
  Matrix3d dR;
  
  //  SO3::Instance().exp(dR, v.head<3>());
  SO3::Instance().cay(dR, v.head<3>());
  
  xb.R = xa.R*dR;
  xb.bg = xa.bg + v.segment<3>(3);
  xb.ba = xa.ba + v.tail<3>();
}


void ImuManifold::dtau(Matrix9d &M, const Vector9d &v)
{
  M.setIdentity();
  Matrix3d dR;
  SO3::Instance().dexp(dR, v.head<3>());
  M.topLeftCorner<3,3>() = dR;
}

void ImuManifold::Adtau(Matrix9d &M, const Vector9d &v)
{
  M.setIdentity();
  Matrix3d g;
  SO3::Instance().exp(g, v.head<3>());
  Matrix3d A;
  SO3::Instance().Ad(A, g);
  M.topLeftCorner<3,3>() = A;
}
