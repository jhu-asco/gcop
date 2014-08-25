#ifndef GCOP_BODY3DMANIFOLD_H
#define GCOP_BODY3DMANIFOLD_H

#include "manifold.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 9, 1> Vector9d;
  typedef Matrix<double, 12, 1> Vector12d;
  typedef Matrix<double, 12, 12> Matrix12d;
  typedef pair<Matrix3d, Vector9d> Body3dState;
  
  class Body3dManifold : public Manifold<Body3dState, 12> {
    
  public:
    static Body3dManifold& Instance();

    void Lift(Vector12d &v,
              const Body3dState &xa,
              const Body3dState &xb);      

    void Retract(Body3dState &xb, 
                 const Body3dState &xa,
                 const Vector12d &v);

    void dtau(Matrix12d &M, const Vector12d &v);

    void Adtau(Matrix12d &M, const Vector12d &v);

  private:
    Body3dManifold();
  };  
}


#endif
