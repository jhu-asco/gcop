#ifndef GCOP_BACOST_H
#define GCOP_BACOST_H

#include <limits>
#include <iostream>
#include "cost.h"
#include "posegraph2d.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, Dynamic, 3> MatrixX3d;
  typedef Matrix<double, 3, Dynamic> Matrix3Xd;

  class BaCost : public Cost<Matrix3d, Vector3d, 3, 3> {
  public:
    BaCost(double tf, const Posegraph2d &pg);

    double Lp(double t, const Matrix3d &x, const Vector3d &u, const VectorXd &p,
              Vector3d *Lx = 0, Matrix3d *Lxx = 0,
              Vector3d *Lu = 0, Matrix3d *Luu = 0,
              Matrix3d *Lxu = 0, 
              VectorXd *Lp = 0, MatrixXd *Lpp = 0, MatrixX3d *Lpx = 0);

    const Posegraph2d &pg;
  };  
}

#endif
