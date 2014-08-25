#ifndef GCOP_PARTICLE2D_H
#define GCOP_PARTICLE2D_H

#include "system.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  /**
   * A simple point-mass (i.e. particle) in 2d with second-order dynamics.
   * The state is
   * \f$ \bf x = (q, v) \in \mathbb{R}^4\f$ with controls \f$ \bf u \in\mathbb{R}^2 \f$ 
   * correspond to turn rate acceleration and forward acceleration.
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */  
  class Particle2d : public System<Vector4d, 4, 2> {
  public:
    Particle2d();
    
    double Step(Vector4d &xb, double t, const Vector4d &xa, 
                const Vector2d &u, double h, const VectorXd *p = 0,
                Matrix<double, 4, 4> *A = 0, Matrix<double, 4, 2> *B = 0, 
                Matrix<double, 4, Dynamic> *C = 0);
    
    double m;  ///< mass
    double r;  ///< radius    
  }; 
}


#endif
