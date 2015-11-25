#ifndef GCOP_GCAR_H
#define GCOP_GCAR_H

#include "system.h"
#include "gcarmanifold.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 4, 1> Vector4d;
  typedef Matrix<double, 4, 4> Matrix4d;
  typedef Matrix<double, 4, 2> Matrix42d;
  typedef Matrix<double, 4, Dynamic> Matrix4Xd;
  
   /**
   * A "geometric" car model with simple dynamics based on the bicycle or Ackerman 
   * steering model. The control system is implemented as an evolution on the Euclidean
   * group SE(2) as opposed to in (x,y,theta) coordinates to avoid issues with orientation
   * algebra, and to use more accurate forward integration. 
   *
   * The state is
   * \f$ \bm x = (g, v) \f$ where \f$ g\in SE(2)\f$ is the pose expressed as a 3x3 matrix, 
   * \f$ v \f$ is the forward velocity, and 
   * the controls \f$ \bm u = (\tau,\tan\phi) \in\mathbb{R}^2\f$ 
   * correspond to rear wheel torque \f$\tau\f$ and to the tan of the steering angle \f$\phi\f$
   *
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  class Gcar : public System<GcarState, 4, 2> {
  public:
    Gcar();
    
    double Step(GcarState &xb, double t, const GcarState &xa, 
                const Vector2d &u, double h, const VectorXd *p,
                Matrix4d *A = 0, Matrix42d *B = 0, Matrix4Xd *C = 0);
    
    double l; ///< distance between axles    
    double r; ///< Gain relating wheel torque to forward acceleration, i.e. \f$ \dot v = r\tau\f$
  };  
}


#endif
