#include "point3dcontroller.h"

using namespace gcop;
using namespace std;
using namespace Eigen;

Point3dController::Point3dController(const Point3d &sys, 
                                     Point3dState *xd, 
                                     Vector3d *ad) : 
  Controller(),
  sys(sys), xd(xd), ad(ad)
{
  Kp.setOnes();
  Kd.setOnes();
}

Point3dController::~Point3dController()
{

}

bool Point3dController::Set(Vector3d &u, double t, const Point3dState &x)
{
  // error in base body configuration
  Vector6d ge;

  // full error in configuration
  Vector3d e = (xd ? x.q - xd->q : x.q);
  Vector3d de = (xd ? x.v - xd->v : x.v);
  
  u = -Kp.cwiseProduct(e) - Kd.cwiseProduct(de);
  // add desired acceleration if provided
  if (ad)
    u = u + *ad;

  return true;
}

