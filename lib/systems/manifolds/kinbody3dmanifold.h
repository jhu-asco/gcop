#ifndef GCOP_KINBODY3DMANIFOLD_H
#define GCOP_KINBODY3DMANIFOLD_H

#include "manifold.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  typedef Matrix<double, 6, 1> Vector6d;
  
  class Kinbody3dManifold : public Manifold<Matrix4d, 6> {
  public:
    
    static Kinbody3dManifold& Instance();
    
    void Lift(Vector6d &v,
              const Matrix4d &xa,
              const Matrix4d &xb);      
    
    void Retract(Matrix4d &xb,
                 const Matrix4d &xa,
                 const Vector6d &v);    

    void dtau(Matrix4d &M, const Vector6d &v);

    void Adtau(Matrix4d &M, const Vector6d &v);
    
  private:
    Kinbody3dManifold();
  };
}

#endif
