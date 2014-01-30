#ifndef GCOP_CAR_H
#define GCOP_CAR_H

#include "system.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 5, 1> Vector5d;
  typedef Matrix<double, 5, 5> Matrix5d;
  typedef Matrix<double, 5, 2> Matrix52d;
  typedef Matrix<double, 5, Dynamic> Matrix5Xd;

   /**
   * A simple rear-drive car model with 2nd order dynamics. 
   *
   * The state is
   * \f$ \bm x = (x,y,\theta, v, k) \in\mathbb{R}^5\f$ and controls are \f$ \bm u = (u_1,u_2) \in\mathbb{R}^5\f$ 
   * correspond to forward acceleration and curvature rate
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  class Car : public System<Vector5d, Vector2d, 5, 2>
  {
  public:
    Car();
    
    double Step(Vector5d &xb, double t, const Vector5d &xa,
                const Vector2d &u, double h, const VectorXd *p = 0,
                Matrix5d *A = 0, Matrix52d *B = 0, Matrix5Xd *C = 0);
    
    double l; ///< distance between axles
  };  
}


#endif
