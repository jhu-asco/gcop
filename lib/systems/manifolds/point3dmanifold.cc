#include <limits>
#include "point3dmanifold.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

Point3dManifold::Point3dManifold() : Manifold()
{
}


Point3dManifold& Point3dManifold::Instance() 
{
  static Point3dManifold instance;
  return instance;
}

void Point3dManifold::Lift(Vector6d &dx, 
                       const Point3dState &xa,
                       const Point3dState &xb) 
{
  dx.head<3>() = xb.q - xa.q;
  dx.tail<3>() = xb.v - xa.v;
}
 

void Point3dManifold::Retract(Point3dState &xb, 
                              const Point3dState &xa,
                              const Vector6d &dx) 
{
  xb.q = xa.q + dx.head<3>();
  xb.v = xa.v + dx.tail<3>();
}
