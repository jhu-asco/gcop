#ifndef GCOP_RCCAR1_H
#define GCOP_RCCAR1_H

#include "system.h"
#include "rccar.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 4, 2> Matrix42d;
  typedef Matrix<double, 4, Dynamic> Matrix4pd;
  //typedef Matrix<double, 6, 1> Vector6d;

   /**
   * A simple rear-drive Rccar model with 2nd order dynamics. 
   *
   * The state is
   * \f$ \bm x = (x,y,\theta, v) \in\mathbb{R}^4\f$ and controls are \f$ \bm u = (a, tan(phi) ) \in\mathbb{R}^5\f$ 
   * correspond to forward acceleration a and steering angle phi
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Rccar1 : public Rccar
  {
  public:
    Rccar1();
    
    double Step(Vector4d &xb, double t, const Vector4d &xa,
                const Vector2d &u, double h, const VectorXd *p,
                Matrix4d *A = 0, Matrix42d *B = 0, Matrix4pd *C = 0);
    inline double bind(const double &input, int index)
    {
      return (input > U.ub[index]? U.ub[index]:(input < U.lb[index])?U.lb[index]:input);
    }

    double msteer, csteer;///<Transform for converting input into steering angle
    double mdrive, cdrive;///<Transform for converting input into desired velocity
    double ktorque;///<Proportional constant for how much torque is applied for matching the current velocity with the desired velocity
  };
}


#endif
