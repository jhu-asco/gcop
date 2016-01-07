#ifndef GCOP_GCARMANIFOLD_H
#define GCOP_GCARMANIFOLD_H

#include "manifold.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 4, 1> Vector4d;
  typedef Matrix<double, 4, 4> Matrix4d;

  class GcarState {
  public:
    GcarState(bool clear = true) {
      if (clear)
        Clear();
    }

    GcarState(Matrix3d gin,double vin):g(gin),v(vin) {
    }

    virtual ~GcarState() {
    }

    void Clear() {
      g.setIdentity();
      v = 0;
    }

    Matrix3d g; ///< SE(2) pose
    double v;   ///< forward velocity   
  };
  
  class GcarManifold : public Manifold<GcarState, 4> {
  public:
    
    static GcarManifold& Instance();
    
    void Lift(Vector4d &dx,
              const GcarState &xa,
              const GcarState &xb);      

    void Retract(GcarState &xb,
                 const GcarState &xa,
                 const Vector4d &dx);    


    void dtau(Matrix4d &M, const Vector4d &v);

    void dtauinv(Matrix4d &M, const Vector4d &v);

    void Adtau(Matrix4d &M, const Vector4d &v);

  private:
    GcarManifold();
  };
}

#endif
