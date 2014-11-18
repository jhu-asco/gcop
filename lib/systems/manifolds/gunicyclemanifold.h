#ifndef GCOP_GUNICYCLEMANIFOLD_H
#define GCOP_GUNICYCLEMANIFOLD_H

#include "manifold.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 5, 1> Vector5d;
  typedef Matrix<double, 5, 5> Matrix5d;
  typedef pair<Matrix3d, Vector2d> M3V2d;
  
  class GunicycleManifold : public Manifold<M3V2d, 5> {
  public:
    
    static GunicycleManifold& Instance();
    
    void Lift(Vector5d &dx,
              const M3V2d &xa,
              const M3V2d &xb);      

    void Retract(M3V2d &xb,
                 const M3V2d &xa,
                 const Vector5d &dx);    

    void dtau(Matrix5d &M, const Vector5d &v);

    void Adtau(Matrix5d &M, const Vector5d &v);

  private:
    GunicycleManifold();
  };
}

#endif
