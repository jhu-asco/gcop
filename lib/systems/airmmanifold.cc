#include <limits>
#include "airmmanifold.h"
#include "so3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

AirmManifold::AirmManifold()
{
}


AirmManifold& AirmManifold::Instance() 
{
  static AirmManifold instance;
  return instance;
}

void AirmManifold::Lift(Vector16d &v, 
                        const AirmState &xa,
                        const AirmState &xb) 
{
  const Matrix3d &Ra = xa.first;  
  const Matrix3d &Rb = xb.first;
  
  Vector3d eR;
  //  SO3::Instance().log(eR, Ra.transpose()*Rb);
  SO3::Instance().cayinv(eR, Ra.transpose()*Rb);
  
  v.head<3>() = eR;
  v.tail<13>() = xb.second - xa.second;
}
 

void AirmManifold::Retract(AirmState &xb, 
                           const AirmState &xa,
                           const Vector16d &v) 
{
  Matrix3d dR;
  
  SO3::Instance().cay(dR, v.head<3>());
  
  xb.first = xa.first*dR;
  xb.second = xa.second + v.tail<13>();
}

