#include <limits>
#include "gunicyclemanifold.h"
#include "se2.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

GunicycleManifold::GunicycleManifold()
{
}


GunicycleManifold& GunicycleManifold::Instance() 
{
  static GunicycleManifold instance;
  return instance;
}


void GunicycleManifold::Lift(Vector5d &v, 
                             const M3V2d &xa,
                             const M3V2d &xb) 
{
  Matrix3d dg;
  SE2::Instance().diff(dg, xa.first, xb.first);
  Vector3d ge;
  SE2::Instance().g2q(ge, dg);
  v.head(3) = ge;
  v.tail(2) = xb.second - xa.second;
}

void GunicycleManifold::Retract(M3V2d &xb, 
                                const M3V2d &xa,
                                const Vector5d &v) 
{
  Matrix3d dg;
  SE2::Instance().q2g(dg, v.head(3));
  xb.first = xa.first*dg;
  xb.second = xa.second + v.tail(2);
}

void GunicycleManifold::dtau(Matrix5d &M, const Vector5d &v)
{
  M.setIdentity();
  double c = cos(v[0]);
  double s = sin(v[0]);
  M(1,1) = c; M(1,2) = -s;
  M(2,1) = s; M(2,2) = c;
}

void GunicycleManifold::Adtau(Matrix5d &M, const Vector5d &v)
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
