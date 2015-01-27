#ifndef GCOP_KINBODY3D_H
#define GCOP_KINBODY3D_H

#include "system.h"
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
 
   /**
   * A simple kinematic control system modeled as a rigid body in 2d. It is meant
   * for parameter estimation problems since there is no dynamics. 
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Kinbody3d : public System<Matrix4d, 6, 6> {
  public:

    Kinbody3d();
    
    double Step(Matrix4d &xb, double t, const Matrix4d &xa,
                const Vector6d &u, double h,  const VectorXd *p = 0,
                Matrix6d *A = 0, Matrix6d *B = 0, Matrix<double, 6, Dynamic> *C = 0);
    
    
    Vector3d d; ///< body dimensions
    
  };
}


#endif
