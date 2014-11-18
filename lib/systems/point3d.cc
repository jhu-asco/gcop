#include <limits>
#include "point3d.h"
#include "point3dmanifold.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;


Point3d::Point3d() : System<Point3dState, 6, 3>(Point3dManifold::Instance()) 
{
  sa = .01;
}


Point3d::~Point3d()
{
}

double Point3d::Step(Point3dState &xb, double t, const Point3dState &xa, 
                     const Vector3d &u, double dt, const VectorXd *p,
                     Matrix6d *A, Matrix63d *B, Matrix6Xd *C) {
  
  xb.q = xa.q + dt*xa.v;
  xb.v = xa.v + dt*u;
  
  // needed for filtering
  if (A) {
    A->setIdentity();
    A->topRightCorner<3,3>().diagonal().setConstant(dt);
  }
  return 0;
}

bool Point3d::Noise(Matrix6d &Q, double t, const Point3dState &x, const Vector3d &u, 
                    double dt, const VectorXd *p) {
  
  double sa2 = sa*sa;
  double dt2 = dt*dt;
  double dt3 = dt2*dt;

  Q.setZero();
  Q.topLeftCorner<3,3>().diagonal().setConstant(sa2*dt3/3);
  Q.topRightCorner<3,3>().diagonal().setConstant(sa2*dt2/2);
  Q.bottomLeftCorner<3,3>().diagonal().setConstant(sa2*dt2/2);
  Q.bottomRightCorner<3,3>().diagonal().setConstant(sa2*dt);
  return true;
}
