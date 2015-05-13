#ifndef GCOP_BODY2DMANIFOLD_H
#define GCOP_BODY2DMANIFOLD_H

#include "manifold.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef pair<Matrix3d, Vector3d> Body2dState;
  
  class Body2dManifold : public Manifold<Body2dState, 6> {
  public:
    
    static Body2dManifold& Instance();
    
    void Lift(Vector6d &dx,
              const Body2dState &xa,
              const Body2dState &xb);      

    void Retract(Body2dState &xb,
                 const Body2dState &xa,
                 const Vector6d &dx);    

    void dtau(Matrix6d &M, const Vector6d &v);

    void Adtau(Matrix6d &M, const Vector6d &v);

  private:
    Body2dManifold();
  };
}

#endif
