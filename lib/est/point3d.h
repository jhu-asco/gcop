#ifndef GCOP_POINT3D_H
#define GCOP_POINT3D_H

#include "system.h"
#include "point3dmanifold.h"
#include <limits>
#include <iostream>
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 3> Matrix63d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 6, Dynamic> Matrix6Xd;
  
  /**
   * An POINT3D system
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Point3d : public System<Point3dState, 6, 3> {        
  public:
  
  Point3d();
  
  virtual ~Point3d();
  
  double Step(Point3dState &xb, double t, const Point3dState &xa, 
              const Vector3d &u, double h, const VectorXd *p = 0,
              Matrix6d *A = 0, Matrix63d *B = 0, Matrix6Xd *C = 0);
  
  bool Noise(Matrix6d &Q, double t, const Point3dState &x, const Vector3d &u, 
             double dt, const VectorXd *p = 0);
  
  double sa; ///< white-noise acceleration stdev

  };  
}

#endif
