#ifndef GCOP_POINT3DMANIFOLD_H
#define GCOP_POINT3DMANIFOLD_H

#include "manifold.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  
  struct Point3dState {
  Point3dState() : 
    q(Vector3d::Zero()), 
      v(Vector3d::Zero()),       
      P(Matrix6d::Identity()) {
    }
    
    Vector3d q;  ///< position
    Vector3d v;  ///< position    
    Matrix6d P;   ///< covariance
  };
  
  //  typedef pair<Matrix3d, Vector6d> Point3dState;
  
  class Point3dManifold : public Manifold<Point3dState, 6> {
    
  public:
    static Point3dManifold& Instance();

    void Lift(Vector6d &dx,
              const Point3dState &xa,
              const Point3dState &xb);      

    void Retract(Point3dState &xb, 
                 const Point3dState &xa,
                 const Vector6d &dx);

  private:
    Point3dManifold();
  };  
}


#endif
