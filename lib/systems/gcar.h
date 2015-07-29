#ifndef GCOP_GCAR_H
#define GCOP_GCAR_H

#include "system.h"
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 4, 1> Vector4d;
  typedef Matrix<double, 4, 4> Matrix4d;
  typedef Matrix<double, 4, 2> Matrix42d;
  typedef Matrix<double, 4, Dynamic> Matrix4Xd;
  
  typedef pair<Matrix3d, double> M3V1d;

   /**
   * A "geometric" car model with 2nd order dynamics. This is often used to
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
  class Gcar : public System<M3V1d, 4, 2> {
  public:
    Gcar();
    
    double Step(M3V1d &xb, double t, const M3V1d &xa, 
                const Vector2d &u, double h, const VectorXd *p,
                Matrix4d *A = 0, Matrix42d *B = 0, Matrix4Xd *C = 0);
    
    double l; ///< distance between axles    
    double r; ///< Gain on the car acc vs wheel torque
  };  
}


#endif
