#ifndef GCOP_KINBODYPROJTRACKCOST_H
#define GCOP_KINBODYPROJTRACKCOST_H

#include <limits>
#include <iostream>
#include "cost.h"
#include "kinbodyprojtrack.h"
#include "kinbody3d.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, Dynamic, 6> MatrixX6d;


  class KinbodyProjTrackCost : public Cost<Matrix4d, 6, 6> {

  private:
    MatrixXd d_xxt(Vector3d x);
    MatrixXd kron(MatrixXd A, MatrixXd B);
    double obsCost(double t, double h, const Matrix4d &x, const VectorXd *p);
    void pokeX(Matrix4d &xb, const Matrix4d &xa, const int i, const double eps);
    void pokeX2(Matrix4d &xb, const Matrix4d &xa, const int i, const int j, const double epsi, const double epsj);
    Matrix6d fdLxx(double t, double h, const Matrix4d &x,  const VectorXd *p);
    Matrix6d fdLx(double t, double h, const Matrix4d &x,  const VectorXd *p);

  public:
    KinbodyProjTrackCost(double tf, const KinbodyProjTrack &pg);
    
    double L(double t, const Matrix4d &x, const Vector6d &u, double h,
             const VectorXd *p,
             Vector6d *Lx = 0, Matrix6d *Lxx = 0,
             Vector6d *Lu = 0, Matrix6d *Luu = 0,
             Matrix<double, 6, 6> *Lxu = 0, 
             VectorXd *Lp = 0, MatrixXd *Lpp = 0, MatrixX6d *Lpx = 0);

    void test_d_xxt();
    void test_kron();
    void test_grads();

    const KinbodyProjTrack &pg;
    bool useFdHess;
    bool useFdGrad;
  };  
}

#endif
