#ifndef GCOP_UNICYCLE_H
#define GCOP_UNICYCLE_H

#include "system.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 5, 1> Vector5d;
  typedef Matrix<double, 5, 5> Matrix5d;
  typedef Matrix<double, 5, 2> Matrix52d;

   /**
   * A standard unicycle model with 2nd order dynamics. This is often used to
   * model simple wheeled robots or fixed-wing aircraft. 
   *
   * The state is
   * \f$ \bm x = (x,y,\theta,v,\omega) \in\mathbb{R}^5\f$ and controls are \f$ \bm u = (u_1,u_2) \in\mathbb{R}^5\f$ 
   * correspond to forward acceleration and turn rate acceleration.
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  class Unicycle : public System<Vector5d, Vector2d, 5, 2>
  {
  public:
    Unicycle();
    
    double Step(Vector5d &xb, double t, const Vector5d &xa, 
             const Vector2d &u, double h,
             Matrix5d *A = 0, Matrix52d *B = 0);
    
    double dx; ///< distance between left and right tires
    double dy; ///< distance between front and back axles
  };  
}


#endif
