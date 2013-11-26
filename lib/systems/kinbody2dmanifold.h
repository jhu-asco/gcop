#ifndef GCOP_KINBODY2DMANIFOLD_H
#define GCOP_KINBODY2DMANIFOLD_H

#include "manifold.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  class Kinbody2dManifold : public Manifold<Matrix3d, 3> {
  public:
    
    static Kinbody2dManifold& Instance();
    
    void Lift(Vector3d &v,
              const Matrix3d &xa,
              const Matrix3d &xb);      
    
    void Retract(Matrix3d &xb,
                 const Matrix3d &xa,
                 const Vector3d &v);    

    void dtau(Matrix3d &M, const Vector3d &v);

    void Adtau(Matrix3d &M, const Vector3d &v);
    
  private:
    Kinbody2dManifold();
  };
}

#endif
