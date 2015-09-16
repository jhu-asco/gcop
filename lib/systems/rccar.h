#ifndef GCOP_RCCAR_H
#define GCOP_RCCAR_H

#include "system.h"
#include <limits>
//#define RCCAR_RADIO_INPUTS

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 4, 2> Matrix42d;
  typedef Matrix<double, 4, Dynamic> Matrix4pd;
  typedef System<Vector4d, 4, 2> BaseSystem;
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
  class Rccar : public BaseSystem 
  {
  public:
    Rccar(int np = 2);
    
    virtual double Step(Vector4d &xb, double t, const Vector4d &xa,
                        const Vector2d &u, double h, const VectorXd *p,
                        Matrix4d *A = 0, Matrix42d *B = 0, Matrix4pd *C = 0);

    /*virtual double Step(Vector4d &xb, double t, const Vector4d &xa,
                        const Vector2d &u, double h, const VectorXd *p,
                        const Vector4d &w,                         
                        Matrix4d *A = 0, Matrix42d *B = 0, Matrix4pd *C = 0, Matrix4d *D = 0);
    */

    virtual void StateAndControlsToFlatAndDerivatives(vector<VectorXd> &y, const Vector4d &x, 
      const std::vector<Vector2d> &u);

		virtual void StateAndControlsToFlat(VectorXd &y, const Vector4d &x, const Vector2d &u);

		virtual void FlatToStateAndControls(Vector4d &x, std::vector<Vector2d> &u, 
      const std::vector<VectorXd> &y);

    double l; ///< distance between axles    
    double r; ///< Gain on the car acc vs wheel torque
    double h; ///< Internal step size
  };
}


#endif
