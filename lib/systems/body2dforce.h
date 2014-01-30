#ifndef GCOP_BODY2DFORCE_H
#define GCOP_BODY2DFORCE_H

#include "force.h"
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef pair<Matrix3d, Vector3d> M3V3d;
  typedef Matrix<double, 3, 6> Matrix36d;

  class Body2dForce : public Force<M3V3d, 3, 6, 3> {
  public:
    Body2dForce(bool fextParam = false);
    
    virtual void Set(Vector3d &f, const M3V3d &x, double t, 
                     const Vector3d &u, double h, const VectorXd *p = 0,
                     Matrix36d *A = 0, Matrix3d *B = 0, Matrix<double, 3, Dynamic> *C = 0);

    Vector3d D;      ///< linear damping terms
    Vector3d fext;   ///< constant external force
    bool fextParam;  ///< treat the (x,y) forces in fext as a parameter
  };
}


#endif
