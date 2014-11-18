#ifndef GCOP_AIRMMANIFOLD_H
#define GCOP_AIRMMANIFOLD_H

#include "manifold.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 13, 1> Vector9d;
  typedef Matrix<double, 16, 1> Vector16d;
  typedef pair<Matrix3d, Vector13d> AirmState;
  
  class AirmManifold : public Manifold<AirmState, 16> {
    
  public:
    
    static AirmManifold& Instance();

    void Lift(Vector16d &v,
              const AirmState &xa,
              const AirmState &xb);      

    void Retract(AirmState &xb, 
                 const AirmState &xa,
                 const Vector16d &v);

  private:
    AirmManifold();
  };  
}


#endif
