#ifndef GCOP_BODY2D_H
#define GCOP_BODY2D_H

#include "system.h"
#include <utility>
#include "body2dforce.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 6, 3> Matrix63d;
  typedef Matrix<double, 3, 6> Matrix36d;
  
  typedef pair<Matrix3d, Vector3d> M3V3d;

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
  class Body2d : public System<M3V3d, Vector3d, 6, 3> {
  public:
    Body2d(Body2dForce *f = 0);

    double Step(M3V3d &xb, double t, const M3V3d &xa,
                const Vector3d &u, double h, const VectorXd *p = 0,
                Matrix6d *A = 0, Matrix63d *B = 0, Matrix<double, 6, Dynamic> *C = 0);
   
    Vector2d d; ///< body dimensions
    Vector3d I; ///< inertia components
    //    Vector3d D; ///< linear damping terms

    Body2dForce *force;  ///< force from controls and external disturbances
    //    Vector3d f; ///< constant external force
  };
}


#endif
