#include <limits>
#include "gcarmanifold.h"
#include "se2.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

GcarManifold::GcarManifold()
{
}


GcarManifold& GcarManifold::Instance() 
{
  static GcarManifold instance;
  return instance;
}


void GcarManifold::Lift(Vector4d &v, 
                        const M3V1d &xa,
                        const M3V1d &xb) 
{
  Matrix3d dg;
  SE2::Instance().diff(dg, xa.first, xb.first);
  Vector3d ge;
  SE2::Instance().g2q(ge, dg);
  v.head(3) = ge;
  v(3) = xb.second - xa.second;
}

void GcarManifold::Retract(M3V1d &xb, 
                           const M3V1d &xa,
                           const Vector4d &v) 
{
  Matrix3d dg;
  SE2::Instance().q2g(dg, v.head(3));
  xb.first = xa.first*dg;
  xb.second = xa.second + v(3);
}

void GcarManifold::dtau(Matrix4d &M, const Vector4d &v)
{
  M.setIdentity();
  double c = cos(v[0]);
  double s = sin(v[0]);
  M(1,1) = c; M(1,2) = -s;
  M(2,1) = s; M(2,2) = c;
}

void GcarManifold::Adtau(Matrix4d &M, const Vector4d &v)
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
