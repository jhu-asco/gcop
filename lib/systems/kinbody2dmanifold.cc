#include "kinbody2dmanifold.h"
#include "se2.h"

using namespace gcop;
using namespace Eigen;

Kinbody2dManifold::Kinbody2dManifold() : Manifold(3)
{
}


Kinbody2dManifold& Kinbody2dManifold::Instance() 
{
  static Kinbody2dManifold instance;
  return instance;
}


void Kinbody2dManifold::Lift(Vector3d &v,
                      const Matrix3d &xa,
                      const Matrix3d &xb)
{
  Matrix3d dg;
  SE2::Instance().diff(dg, xa, xb);
  SE2::Instance().g2q(v, dg);
}


void Kinbody2dManifold::Retract(Matrix3d &xb,
                                const Matrix3d &xa,
                                const Vector3d &v)
{
  Matrix3d dg;
  SE2::Instance().q2g(dg, v);
  xb = xa*dg;
}


void Kinbody2dManifold::dtau(Matrix3d &M, const Vector3d &v)
{
  M.setIdentity();
  double c = cos(v[0]);
  double s = sin(v[0]);
  M(1,1) = c; M(1,2) = -s;
  M(2,1) = s; M(2,2) = c;
}


void Kinbody2dManifold::Adtau(Matrix3d &M, const Vector3d &v)
{
  const double c = cos(v[0]);
  const double s = sin(v[0]);
  const double &x = v[1];
  const double &y = v[2];

  M.setIdentity();
  M(1,0) = y;
  M(1,1) = c;
  M(1,2) = -s;
  M(2,0) = -x;
  M(2,1) = s;
  M(2,2) = c;
}
