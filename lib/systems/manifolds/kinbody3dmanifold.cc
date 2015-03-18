#include "kinbody3dmanifold.h"
#include "se3.h"

using namespace gcop;
using namespace Eigen;

Kinbody3dManifold::Kinbody3dManifold() : Manifold(6)
{
}


Kinbody3dManifold& Kinbody3dManifold::Instance() 
{
  static Kinbody3dManifold instance;
  return instance;
}


void Kinbody3dManifold::Lift(Vector6d &v,
                      const Matrix4d &xa,
                      const Matrix4d &xb)
{
  Matrix4d dg;
  SE3::Instance().diff(dg, xa, xb);
  SE3::Instance().g2q(v, dg);
}


void Kinbody3dManifold::Retract(Matrix4d &xb,
                                const Matrix4d &xa,
                                const Vector6d &v)
{
  Matrix4d dg;
  SE3::Instance().q2g(dg, v);
  xb = xa*dg;
}


void Kinbody3dManifold::dtau(Matrix4d &M, const Vector6d &v)
{
  M.setIdentity();
  Matrix3d dR;
  SO3::Instance().dcay(dR, v.head<3>());
  M.topLeftCorner<3,3>() = dR;
}


void Kinbody3dManifold::Adtau(Matrix4d &M, const Vector6d &v)
{
  M.setIdentity();
  Matrix3d g;
  SO3::Instance().cay(g, v.head<3>());
  Matrix3d A;
  SO3::Instance().Ad(A, g);
  M.topLeftCorner<3,3>() = A;
}
