#ifndef GCOP_BODY3DTRACKCOST_H
#define GCOP_BODY3DTRACKCOST_H

#include <limits>
#include <iostream>
#include "cost.h"
#include "body3dtrack.h"
#include "body3d.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, Dynamic, 12> MatrixX12d;
  typedef Matrix<double, 12, Dynamic> Matrix12Xd;
  
  class Body3dTrackCost : public Cost<Body3dState, 12, 6> {
  public:
    Body3dTrackCost(double tf, const Body3dTrack &pg);
    
    double L(double t, const Body3dState &x, const Vector6d &u, double h,
             const VectorXd *p,
             Vector12d *Lx = 0, Matrix12d *Lxx = 0,
             Vector6d *Lu = 0, Matrix6d *Luu = 0,
             Matrix<double, 12, 6> *Lxu = 0, 
             VectorXd *Lp = 0, MatrixXd *Lpp = 0, MatrixX12d *Lpx = 0);

    const Body3dTrack &pg;
  };  
}

#endif
