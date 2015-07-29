#ifndef GCOP_GCARMANIFOLD_H
#define GCOP_GCARMANIFOLD_H

#include "manifold.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 4, 1> Vector4d;
  typedef Matrix<double, 4, 4> Matrix4d;
  typedef pair<Matrix3d, double> M3V1d;
  
  class GcarManifold : public Manifold<M3V1d, 4> {
  public:
    
    static GcarManifold& Instance();
    
    void Lift(Vector4d &dx,
              const M3V1d &xa,
              const M3V1d &xb);      

    void Retract(M3V1d &xb,
                 const M3V1d &xa,
                 const Vector4d &dx);    

    void dtau(Matrix4d &M, const Vector4d &v);

    void Adtau(Matrix4d &M, const Vector4d &v);

  private:
    GcarManifold();
  };
}

#endif
