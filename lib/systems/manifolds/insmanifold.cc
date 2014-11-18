#include <limits>
#include "insmanifold.h"
#include "so3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

InsManifold::InsManifold() : Manifold()
{
}


InsManifold& InsManifold::Instance() 
{
  static InsManifold instance;
  return instance;
}

void InsManifold::Lift(Vector15d &v, 
                       const InsState &xa,
                       const InsState &xb) 
{
  Vector3d eR;
  //  SO3::Instance().log(eR, xa.R.transpose()*xb.R);
  SO3::Instance().cayinv(eR, xa.R.transpose()*xb.R);
  
  v.head<3>() = eR;
  v.segment<3>(3) = xb.bg - xa.bg;
  v.segment<3>(6) = xb.ba - xa.ba;
  v.segment<3>(9) = xb.p - xa.p;
  v.segment<3>(12) = xb.v - xa.v;
}
 

void InsManifold::Retract(InsState &xb, 
                          const InsState &xa,
                          const Vector15d &v) 
{
  Matrix3d dR;
  
  //  SO3::Instance().exp(dR, v.head<3>());
  SO3::Instance().cay(dR, v.head<3>());
  
  xb.R = xa.R*dR;
  xb.bg = xa.bg + v.segment<3>(3);
  xb.ba = xa.ba + v.segment<3>(6);
  xb.p = xa.p + v.segment<3>(9);
  xb.v = xa.v + v.segment<3>(12);
}


void InsManifold::dtau(Matrix15d &M, const Vector15d &v)
{
  M.setIdentity();
  Matrix3d dR;
  SO3::Instance().dexp(dR, v.head<3>());
  M.topLeftCorner<3,3>() = dR;
}

void InsManifold::Adtau(Matrix15d &M, const Vector15d &v)
{
  M.setIdentity();
  Matrix3d g;
  SO3::Instance().exp(g, v.head<3>());
  Matrix3d A;
  SO3::Instance().Ad(A, g);
  M.topLeftCorner<3,3>() = A;
}
