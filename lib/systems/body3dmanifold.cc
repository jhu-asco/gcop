#include <limits>
#include "body3dmanifold.h"
#include "so3.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Body3dManifold::Body3dManifold()
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
  const Matrix3d &Ra = xa.first;  
  const Matrix3d &Rb = xb.first;
  
  Vector3d eR;
  //  SO3::Instance().log(eR, Ra.transpose()*Rb);
  SO3::Instance().cayinv(eR, Ra.transpose()*Rb);
  
  v.head<3>() = eR;
  v.tail<9>() = xb.second - xa.second;
}
 

void Body3dManifold::Retract(Body3dState &xb, 
                             const Body3dState &xa,
                             const Vector12d &v) 
{
  Matrix3d dR;
  
  SO3::Instance().cay(dR, v.head<3>());
  
  xb.first = xa.first*dR;
  xb.second = xa.second + v.tail<9>();
}


void Body3dManifold::dtau(Matrix12d &M, const Vector12d &v)
{
  M.setIdentity();
  Matrix3d dR;
  SO3::Instance().dcay(dR, v.head<3>());
  M.topLeftCorner<3,3>() = dR;
}

void Body3dManifold::Adtau(Matrix12d &M, const Vector12d &v)
{
  M.setIdentity();
  Matrix3d g;
  SO3::Instance().cay(g, v.head<3>());
  Matrix3d A;
  SO3::Instance().Ad(A, g);
  M.topLeftCorner<3,3>() = A;
}
