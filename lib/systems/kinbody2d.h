#ifndef GCOP_KINBODY2D_H
#define GCOP_KINBODY2D_H

#include "system.h"
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
   /**
   * A simple kinematic control system modeled as a rigid body in 2d. It is meant
   * for parameter estimation problems since there is no dynamics. 
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Kinbody2d : public System<Matrix3d, 3, 3> {
  public:

    Kinbody2d();
    
    double Step(Matrix3d &xb, double t, const Matrix3d &xa,
                const Vector3d &u, double h,  const VectorXd *p,
                Matrix3d *A = 0, Matrix3d *B = 0, Matrix<double, 3, Dynamic> *C = 0);
    
    
    Vector2d d; ///< body dimensions
    
  };
}


#endif
