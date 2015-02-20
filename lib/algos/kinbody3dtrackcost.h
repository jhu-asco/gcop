#ifndef GCOP_KINBODY3DTRACKCOST_H
#define GCOP_KINBODY3DTRACKCOST_H

#include <limits>
#include <iostream>
#include "cost.h"
#include "kinbody3dtrack.h"
#include "kinbody3d.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, Dynamic, 6> MatrixX6d;


  class Kinbody3dTrackCost : public Cost<Matrix4d, 6, 6> {
  public:
    Kinbody3dTrackCost(double tf, const Kinbody3dTrack &pg);
    
    double L(double t, const Matrix4d &x, const Vector6d &u, double h,
             const VectorXd *p,
             Vector6d *Lx = 0, Matrix6d *Lxx = 0,
             Vector6d *Lu = 0, Matrix6d *Luu = 0,
             Matrix<double, 6, 6> *Lxu = 0, 
             VectorXd *Lp = 0, MatrixXd *Lpp = 0, MatrixX6d *Lpx = 0);

    double obsCost(double t, double h, const Matrix4d &x, const VectorXd *p);
    void pokeX(Matrix4d &xb, const Matrix4d &xa, const int i, const double eps);
    void pokeX2(Matrix4d &xb, const Matrix4d &xa, const int i, const int j, const double epsi, const double epsj);
    Matrix6d fdLxx(double t, double h, const Matrix4d &x,  const VectorXd *p);

    const Kinbody3dTrack &pg;
    bool useFdHess;
  };  
}

#endif
