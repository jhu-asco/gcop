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


  template <int _nu = 6>
  class KinbodyProjTrackCost : public Cost<Matrix4d, 6, _nu> {

  private:
    MatrixXd d_xxt(Vector3d x);
    MatrixXd kron(MatrixXd A, MatrixXd B);
  public:
    typedef Matrix<double, _nu, _nu> Matrixud;
    typedef Matrix<double, _nu, 1> Vectorud;
    typedef Matrix<double, 6, 1> Vector6d;
    typedef Matrix<double, 6, 6> Matrix6d;

    KinbodyProjTrackCost(double tf, const KinbodyProjTrack<_nu> &pg);
    
    double L(double t, const Matrix4d &x, const Vectorud &u, double h,
             const VectorXd *p,
             Vector6d *Lx = 0, Matrix6d *Lxx = 0,
             Vectorud *Lu = 0, Matrixud *Luu = 0,
             Matrix<double, 6, _nu> *Lxu = 0, 
             VectorXd *Lp = 0, MatrixXd *Lpp = 0, MatrixX6d *Lpx = 0);

    void test_d_xxt();
    void test_kron();
    void test_grads();

    const KinbodyProjTrack<_nu> &pg;
  };  

template <int _nu>
KinbodyProjTrackCost<_nu>::KinbodyProjTrackCost(double tf, const KinbodyProjTrack<_nu> &pg) :
  Cost<Matrix4d, 6, _nu>(pg.sys, tf), pg(pg)
{
  
}


static Matrix3d r3hat(const Vector3d &a)
{
  SO3& so3 = SO3::Instance();
  Matrix3d a_hat;

  so3.hat(a_hat, a);
  //a_hat << 0,  -a(2),  a(1),
  //       a(2),    0,  -a(0),
  //      -a(1),  a(0),    0;
  return a_hat;
}


static Vector3d cross3(Vector3d a, Vector3d b)
{
  return a.cross(b);
}

template <int _nu>
void KinbodyProjTrackCost<_nu>::test_kron()
{
  MatrixXd A = MatrixXd::Identity(6,6);
  MatrixXd B = MatrixXd::Zero(2,1);
  B << 3, 2;

  cout << "A:" << endl << A << endl;
  cout << "B:" << endl << B << endl;
  cout << "kron(A,B): " << endl << kron(A,B) << endl;
  cout << "kron(B,A): " << endl << kron(B,A) << endl;
}

template <int _nu>
void KinbodyProjTrackCost<_nu>::test_d_xxt()
{
  Vector3d x(3,-2,1);
  cout << "x: " << endl << x << endl; 
  cout << "d_xxt: " << endl << d_xxt(x) << endl; 
}

template <int _nu>
void KinbodyProjTrackCost<_nu>::test_grads()
{
  Vector6d Lx;
  Matrix6d Lxx;
  VectorXd Lp;
  Lp.resize(15,1);
  MatrixXd  Lpp;
  Lpp.resize(15,15); 
  MatrixX6d Lpx;
  Lpx.resize(15,6);

  Lx.setZero();
  Lxx.setZero();
  Lp.setZero();
  Lpp.setZero();
  Lpx.setZero();

  Matrix<double, 5, 3> ls;
  Matrix<double, 5, 3> zs;

  ls << -8.049191900011810,  -4.430035622659032,   0.937630384099677,
   9.150136708685952,   9.297770703985531,  -6.847738366449034,
   9.411855635212312,   9.143338964858913,  -0.292487025543176,
   6.005609377776004,  -7.162273227455693,  -1.564774347474501,
   8.314710503781342,   5.844146591191087,   9.189848527858061;

  zs << 0.438021624604109,  -0.288302180913201,  -0.851480421888765,
   0.234909466397926,   0.387309344924147,  -0.891520618904055,
   0.399428923259037,   0.397492249355858,  -0.826109222177157,
   0.544326983672147,   0.068898439211757,  -0.836038958374887,
   0.705665941886738,   0.264747508473682,  -0.657224722007686;
  

  double L = 0;

  double t = 5.119059895756811;
  double t2 = 5.691258590395267; 
  Matrix3d R1, R2;
  R1 << cos(t), -sin(t), 0,
        sin(t), cos(t),  0,
             0,      0,  1;
  R2 <<      1,      0,  0,
         0, cos(t2), -sin(t2),
         0, sin(t2),  cos(t2);
  Matrix3d R = R2*R1; // orientation
  Vector3d xp;     // position
  xp <<   -4.920527348259758+5,
  26.535034245560773,
  15.294369849016380;
  
  cout << "R: " << endl << R << endl;
  cout << "xp: " << endl << xp << endl;
  cout << "ls: " << endl << ls << endl;
  cout << "zs: " << endl << zs << endl;
  cout << "sigmaR: " << endl << 1./pg.cp << endl;

  Matrix3d eye = MatrixXd::Identity(3,3);
   
  for (int i = 0; i < ls.rows(); ++i) {
    const Vector3d &z = zs.block<1,3>(i,0);
    const Vector3d &pf = ls.block<1,3>(i,0);
    Vector3d y = R.transpose()*(pf - xp);
    y /= y.norm();
    Vector3d r = (y - z)/pg.cp;

    L += r.dot(y - z)/2;         // feature cost

    double normlp = (pf - xp).norm();
    Vector3d lphat = (pf - xp)/normlp;
    Matrix3d lphat_d = (eye - lphat*lphat.transpose())/normlp;

    //cout << "feature cost " << i << "=" << r.dot(y - z)/2 << endl;

    Lx.segment<3>(0) += - r3hat(y)*r;
    Lx.segment<3>(3) += - R.transpose()*lphat_d*R*r;

    Lxx.block<3,3>(0,0) += (r3hat(y).transpose()*r3hat(y))/pg.cp + (r3hat(r)*r3hat(y))/2. + (r3hat(r)*r3hat(y)).transpose()/2.;
    Lxx.block<3,3>(3,0) += -r3hat(R.transpose()*lphat_d*R*r)/2. + R.transpose()*lphat_d*R*r3hat(r) - R.transpose() *lphat_d*R*r3hat(y)/pg.cp;
    Lxx.block<3,3>(3,3) += -R.transpose()*(-lphat_d*R*R.transpose()*lphat_d*R/pg.cp + (kron((R*r).transpose(), eye)*d_xxt(pf-xp)*R)/normlp + (eye - lphat*lphat.transpose())*(R*r*lphat.transpose()*R)/(normlp*normlp));

    Lp.segment<3>(3*i) = lphat_d*R*r;

    Lpp.block<3,3>(3*i, 3*i) = lphat_d*R*R.transpose()*lphat_d/pg.cp - (kron((R*r).transpose(), eye)*d_xxt(pf-xp))/normlp - (eye - lphat*lphat.transpose())*(R*r*lphat.transpose())/(normlp*normlp);
    
    Lpx.block<3,3>(3*i, 0) = -lphat_d*R*r3hat(r) + lphat_d*R*r3hat(y)/pg.cp;
    Lpx.block<3,3>(3*i, 3) = -lphat_d*R*R.transpose()*lphat_d*R/pg.cp + kron((R*r).transpose(), eye)*d_xxt(pf-xp)*R/normlp + (eye - lphat*lphat.transpose())*(R*r*lphat.transpose()*R)/(normlp*normlp);
  }

  Lxx.block<3,3>(0,3) = Lxx.block<3,3>(3,0).transpose();

  cout << "L: " << endl << L << endl;
  cout << "Lx: " << endl << Lx << endl;
  cout << "Lxx: " << endl << Lxx << endl;
  cout << "Lp: " << endl << Lp << endl;
  cout << "Lpp: " << endl << Lpp << endl;
  cout << "Lpx: " << endl << Lpx << endl;
}

template <int _nu>
MatrixXd KinbodyProjTrackCost<_nu>::kron(MatrixXd A, MatrixXd B)
{
  int m = B.rows();
  int n = B.cols();
  MatrixXd k = Eigen::MatrixXd(A.rows()*m, A.cols()*n);
  for(int i = 0; i < A.rows(); i++)
  {
    for(int j = 0; j < A.cols(); j++)
    {
      k.block(i*m, j*n, m, n) = A(i,j)*B;
    }
  }
  return k;
}

// Gives Hessian of stack(xhat*xhat^T) with respect to x
template <int _nu>
MatrixXd KinbodyProjTrackCost<_nu>::d_xxt(Vector3d x)
{
  double x1 = x(0);
  double x2 = x(1);
  double x3 = x(2);

  MatrixXd hess = MatrixXd(9,3);

  double mag2_x = x1*x1 + x2*x2 + x3*x3; 

  hess <<                  (2*x1)/mag2_x - (2*x1*x1*x1)/pow(mag2_x,2),   
                                          -(2*x1*x1*x2)/pow(mag2_x,2),
                                          -(2*x1*x1*x3)/pow(mag2_x,2),
                               x2/mag2_x - (2*x1*x2*x1)/pow(mag2_x,2),   
                               x1/mag2_x - (2*x1*x2*x2)/pow(mag2_x,2),  
                                          -(2*x1*x2*x3)/pow(mag2_x,2),
                               x3/mag2_x - (2*x1*x3*x1)/pow(mag2_x,2),
                                          -(2*x1*x3*x2)/pow(mag2_x,2),
                               x1/mag2_x - (2*x1*x3*x3)/pow(mag2_x,2),
                               x2/mag2_x - (2*x1*x2*x1)/pow(mag2_x,2),
                               x1/mag2_x - (2*x1*x2*x2)/pow(mag2_x,2),
                                          -(2*x1*x2*x3)/pow(mag2_x,2),
                                          -(2*x2*x2*x1)/pow(mag2_x,2),
                           (2*x2)/mag2_x - (2*x2*x2*x2)/pow(mag2_x,2),
                                          -(2*x2*x2*x3)/pow(mag2_x,2),
                                          -(2*x2*x3*x1)/pow(mag2_x,2),
                               x3/mag2_x - (2*x2*x3*x2)/pow(mag2_x,2),
                               x2/mag2_x - (2*x2*x3*x3)/pow(mag2_x,2),
                               x3/mag2_x - (2*x1*x3*x1)/pow(mag2_x,2),
                                          -(2*x1*x3*x2)/pow(mag2_x,2),
                               x1/mag2_x - (2*x1*x3*x3)/pow(mag2_x,2),
                                          -(2*x2*x3*x1)/pow(mag2_x,2),
                               x3/mag2_x - (2*x2*x3*x2)/pow(mag2_x,2),
                               x2/mag2_x - (2*x2*x3*x3)/pow(mag2_x,2),
                                          -(2*x3*x3*x1)/pow(mag2_x,2),
                                          -(2*x3*x3*x2)/pow(mag2_x,2),
                           (2*x3)/mag2_x - (2*x3*x3*x3)/pow(mag2_x,2);
  return hess;
}

template <int _nu>
double KinbodyProjTrackCost<_nu>::L(double t, const Matrix4d &x, const Vectorud &u, 
                          double h,
                          const VectorXd *p,
                          Vector6d *Lx, Matrix6d *Lxx,
                          Vectorud *Lu, Matrixud *Luu,
                          Matrix<double, 6, _nu> *Lxu,
                          VectorXd *Lp, MatrixXd *Lpp, MatrixX6d *Lpx)
{
  double L = 0;

  //cout << "x: " << endl << x << endl;

  if (Lu)
    Lu->setZero();

  if (Luu)
    Luu->setZero();

  if (Lxu)
    Lxu->setZero();
  
  if (Lx)
    Lx->setZero();

  if (Lxx)
    Lxx->setZero();

  if (Lp)
    Lp->setZero();

  if (Lpp)
    Lpp->setZero();

  if (Lpx)
    Lpx->setZero();
  
  //const Matrix3d &g = x.first;

  const Matrix3d &R = x.block<3,3>(0,0); // orientation
  const Vector3d &xp = x.block<3,1>(0,3);     // position
  
  int N = pg.Is.size() - 1;
  
  h = this->tf/N;

  //  int k = (h < 1e-16 ? N : (int)round(t/h));

  int k = (int)round(t/h);
  assert(k >=0 && k <= N);
  const vector< pair<int, Vector3d> > &I = pg.Is[k];

  //cout << "Kinbody3dtrackCost: k=" << k << " " << I.size() << endl;

  Matrix3d eye = MatrixXd::Identity(3,3);

  int i0 = 3*pg.extforce;
  for (int i = 0; i < I.size(); ++i) {
    int l = pg.cis[I[i].first];  // feature index
    const Vector3d &z = I[i].second;

    const Vector3d &pf = p->segment<3>(i0 + 3*l);
    Vector3d y = R.transpose()*(pf - xp);
    y /= y.norm();
    Vector3d r = (y - z)/pg.cp;

    L += r.dot(y - z)/2;         // feature cost

    double normlp = (pf - xp).norm();
    Vector3d lphat = (pf - xp)/normlp;
    Matrix3d lphat_d = (eye - lphat*lphat.transpose())/normlp;

    //cout << "feature cost " << i << "=" << r.dot(y - z)/2 << endl;

    if (Lx) {
      //Lx->segment<3>(0) = Lx->segment<3>(0) - cross3(y, r);
      Lx->segment<3>(0) = Lx->segment<3>(0) - r3hat(y)*r;
      Lx->segment<3>(3) = Lx->segment<3>(3) - R.transpose()*lphat_d*R*r;
    }

    if (Lxx) {
      Lxx->block<3,3>(0,0) += (r3hat(y).transpose()*r3hat(y))/pg.cp + (r3hat(r)*r3hat(y))/2. + (r3hat(r)*r3hat(y)).transpose()/2.;
      Lxx->block<3,3>(3,0) += -r3hat(R.transpose()*lphat_d*R*r)/2. + R.transpose()*lphat_d*R*r3hat(r) - R.transpose() *lphat_d*R*r3hat(y)/pg.cp;
      Lxx->block<3,3>(3,3) += -R.transpose()*(-lphat_d*R*R.transpose()*lphat_d*R/pg.cp + (kron((R*r).transpose(), eye)*d_xxt(pf-xp)*R)/normlp + (eye - lphat*lphat.transpose())*(R*r*lphat.transpose()*R)/(normlp*normlp));
    }

    if (Lp)
      Lp->segment<3>(i0 + 3*l) = lphat_d*R*r;

    if (Lpp)
    {
      Lpp->block<3,3>(i0 + 3*l, i0 + 3*l) = lphat_d*R*R.transpose()*lphat_d/pg.cp - (kron((R*r).transpose(), eye)*d_xxt(pf-xp))/normlp - (eye - lphat*lphat.transpose())*(R*r*lphat.transpose())/(normlp*normlp);
    }

    if (Lpx) {
      Lpx->block<3,3>(i0 + 3*l, 0) = -lphat_d*R*r3hat(r) + lphat_d*R*r3hat(y)/pg.cp;
      Lpx->block<3,3>(i0 + 3*l, 3) = -lphat_d*R*R.transpose()*lphat_d*R/pg.cp + kron((R*r).transpose(), eye)*d_xxt(pf-xp)*R/normlp + (eye - lphat*lphat.transpose())*(R*r*lphat.transpose()*R)/(normlp*normlp);
      //Lpx->block<3,3>(i0 + 3*l, 3) = Eigen::MatrixXd::Zero(3,3);
    }
  }

  if(Lxx)
  {
    Lxx->block<3,3>(0,3) = Lxx->block<3,3>(3,0).transpose();
  }

  if (pg.forces && k < N) {
    Vectorud du = u - pg.uos[k];
    Vectorud wdu = du.cwiseQuotient(pg.cw);

    L += du.dot(wdu)/2;
    if (Lu)
      *Lu = wdu;
    if (Luu)
      *Luu = Vector6d::Ones().cwiseQuotient(pg.cw).asDiagonal();
  }

  /*
  if(Lp)
    cout << "Lp: " << endl << *Lp << endl;
  if(Lpp)
    cout << "Lpp: " << endl << *Lpp << endl;
  if(Lpx)
    cout << "Lpx: " << endl << *Lpx << endl;
  if(Lx)
    cout << "Lx: " << endl << *Lx << endl;
  if(Lxx)
    cout << "Lxx: " << endl << *Lxx << endl;
  */

  return L;
}
}

#endif
