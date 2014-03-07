#ifndef GCOP_GUNICYCLE_H
#define GCOP_GUNICYCLE_H

#include "system.h"
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 5, 1> Vector5d;
  typedef Matrix<double, 5, 5> Matrix5d;
  typedef Matrix<double, 5, 2> Matrix52d;
  typedef Matrix<double, 5, Dynamic> Matrix5Xd;
  
  typedef pair<Matrix3d, Vector2d> M3V2d;

   /**
   * A "geometric" unicycle model with 2nd order dynamics. This is often used to
   * model simple wheeled robots or fixed-wing aircraft. The control system
   * is implemented as an evolution on SE(2) as opposed to a vector space.
   *
   * The state is
   * \f$ \bm x = (g, \omega, v) \f$ where \f$ g\in SE(2)\f$ is the pose, 
   * the controls \f$ \bm u = (u_\omega,u_v) \in\mathbb{R}^2\f$ 
   * correspond to turn rate acceleration and forward acceleration.
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  class Gunicycle : public System<M3V2d, Vector2d, 5, 2> {
  public:
    Gunicycle();
    
    double Step(M3V2d &xb, double t, const M3V2d &xa, 
                const Vector2d &u, double h, const VectorXd *p,
                Matrix5d *A = 0, Matrix52d *B = 0, Matrix5Xd *C = 0);
    
    double dx; ///< distance between left and right tires
    double dy; ///< distance between front and back axles
  };  
}


#endif
