#ifndef GCOP_GUNICYCLECOST_H
#define GCOP_GUNICYCLECOST_H

#include "lqcost.h"
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 5, 1> Vector5d;
  typedef Matrix<double, 5, 5> Matrix5d;
  typedef Matrix<double, 5, 2> Matrix52d;
  typedef pair<Matrix3d, Vector2d> M3V2d;
  

  class GunicycleCost : public LqCost<M3V2d, 5, 2> {
  public:      
    
    GunicycleCost(double tf, const M3V2d &xf);
    
    double L(double t, const M3V2d &x, const Vector2d &u, double h,
             Vector5d *Lx = 0, Matrix5d *Lxx = 0,
             Vector2d *Lu = 0, Matrix2d *Luu = 0,
             Matrix52d *Lxu = 0);
  };
}

#endif
