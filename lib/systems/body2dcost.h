#ifndef GCOP_BODY2DCOST_H
#define GCOP_BODY2DCOST_H

#include "lqcost.h"
#include <limits>
#include "body2dtrack.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 6, 3> Matrix63d;
  typedef pair<Matrix3d, Vector3d> M3V3d;
  
  
  class Body2dCost : public LqCost<M3V3d, 6, 3> {
  public:      
    
    Body2dCost(double tf, const M3V3d &xf);
    
    double L(double t, const M3V3d &x, const Vector3d &u, double h,
             Vector6d *Lx = 0, Matrix6d *Lxx = 0,
             Vector3d *Lu = 0, Matrix3d *Luu = 0,
             Matrix63d *Lxu = 0);

    
    Body2dTrack *track;
    double ko;  ///< obstacle avoidance gain
    
  };
}

#endif
