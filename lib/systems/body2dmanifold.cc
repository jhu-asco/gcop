#include <limits>
#include "body2dmanifold.h"
#include "se2.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Body2dManifold::Body2dManifold()
{
}


Body2dManifold& Body2dManifold::Instance() 
{
  static Body2dManifold instance;
  return instance;
}


void Body2dManifold::Lift(Vector6d &v, 
                          const M3V3d &xa,
                          const M3V3d &xb) 
{
  Matrix3d dg;
  SE2::Instance().diff(dg, xa.first, xb.first);
  Vector3d ge;
  SE2::Instance().g2q(ge, dg);
  v.head(3) = ge;
  v.tail(3) = xb.second - xa.second;
}

void Body2dManifold::Retract(M3V3d &xb, 
                             const M3V3d &xa,
                             const Vector6d &v) 
{
  Matrix3d dg;
  SE2::Instance().q2g(dg, v.head(3));
  xb.first = xa.first*dg;
  xb.second = xa.second + v.tail(3);
}

void Body2dManifold::dtau(Matrix6d &M, const Vector6d &v)
{
  M.setIdentity();
  double c = cos(v[0]);
  double s = sin(v[0]);
  M(1,1) = c; M(1,2) = -s;
  M(2,1) = s; M(2,2) = c;
}

void Body2dManifold::Adtau(Matrix6d &M, const Vector6d &v)
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
