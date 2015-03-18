#ifndef GCOP_KINRCCAR_H
#define GCOP_KINRCCAR_H

#include "kinbody3d.h"
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
  typedef Matrix<double, 6, 2> Matrix62d;
  typedef Matrix<double, 6, 6> Matrix6d;

   /**
   * A simple kinematic car model in SE(3). 
   *
   * Author: Matt Sheckells mshecke1(at)jhu.edu
   */
  class KinRccar : public Kinbody3d<2> {
  public:

    KinRccar();
    
    double Step(Matrix4d &xb, double t, const Matrix4d &xa,
                const Vector2d &u, double h,  const VectorXd *p = 0,
                Matrix6d *A = 0, Matrix62d *B = 0, Matrix<double, 6, Dynamic> *C = 0);
    
    
  };
}


#endif
