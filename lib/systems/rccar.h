#ifndef GCOP_RCCAR_H
#define GCOP_RCCAR_H

#include "system.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 4, 2> Matrix42d;
  typedef Matrix<double, 4, Dynamic> Matrix4Xd;

   /**
   * A simple rear-drive Rccar model with 2nd order dynamics. 
   *
   * The state is
   * \f$ \bm x = (x,y,\theta, v) \in\mathbb{R}^4\f$ and controls are \f$ \bm u = (a, tan(phi) ) \in\mathbb{R}^5\f$ 
   * correspond to forward acceleration a and steering angle phi
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Rccar : public System<Vector4d, Vector2d, 4, 2>
  {
  public:
    Rccar();
    
    double Step(Vector4d &xb, double t, const Vector4d &xa,
                const Vector2d &u, double h, const VectorXd *p = 0,
                Matrix4d *A = 0, Matrix42d *B = 0, Matrix4Xd *C = 0);

    double l; ///< distance between axles    
  };
}


#endif
