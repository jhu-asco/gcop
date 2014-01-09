#include "se3.h"
#include <assert.h>

using namespace gcop;
using namespace Eigen;
using namespace std;

//vec SE3::e[3] = {vec("1,0,0"), vec("0,1,0"),vec("0,0,1")};
//mat SE3::E[3] = {SE3::ad(e[0]), SE3::ad(e[1]), SE3::ad(e[2])};

SE3::SE3() : Group(),
             tol(1e-16),
             so3(SO3::Instance())
{
}

SE3& SE3::Instance()
{
  static SE3 se3;
  return se3;
}



void SE3::inv(Matrix4d &gi, const Matrix4d &g) const
{
  const Matrix3d &R = g.topLeftCorner<3,3>();
  const Vector3d &p = g.topRightCorner<3,1>();

  gi.topLeftCorner<3,3>() = R.transpose();
  gi.topRightCorner<3,1>() = -R.transpose()*p;
  gi.bottomLeftCorner<1,3>().setZero();
  gi(3,3) = 1;
}


void SE3::Ad(Matrix6d &M, const Matrix4d &g) const
{
  const Matrix3d &R = g.topLeftCorner<3,3>();
  const Vector3d &p = g.topRightCorner<3,1>();

  Matrix3d adp;  
  so3.ad(adp, p);

  M.topLeftCorner<3,3>() = R;
  M.topRightCorner<3,3>().setZero();
  M.bottomLeftCorner<3,3>() = adp*R;
  M.bottomRightCorner<3,3>() = R;
}


void SE3::hat(Matrix4d &vh, const Vector6d &v) const
{
  Matrix3d wh;  
  so3.ad(wh, v.head<3>());

  vh.topLeftCorner<3,3>() = wh;
  vh.topRightCorner<3,1>() = v.tail<3>();  
  vh.row(3).setZero();
}


void SE3::hatinv(Vector6d& v, const Matrix4d& vh) const
{
  Vector3d w;
  so3.hatinv(w, vh.topLeftCorner<3,3>());

  v.head<3>() = w;
  v.tail<3>() = vh.topRightCorner<3,1>();
}

void SE3::Tg(Matrix6d &M, const Matrix4d &g) const
{
  const Matrix3d &R = g.topLeftCorner<3,3>();

  M.topLeftCorner<3,3>() = R;
  M.topRightCorner<3,3>().setZero();
  M.bottomLeftCorner<3,3>().setZero();
  M.bottomRightCorner<3,3>() = R;
}


void SE3::ad(Matrix6d &M, const Vector6d &v) const
{
  Matrix3d wh;  
  so3.ad(wh, v.head<3>());

  Matrix3d vh;  
  so3.ad(vh, v.tail<3>());

  M.topLeftCorner<3,3>() = wh;
  M.topRightCorner<3,3>().setZero();
  M.bottomLeftCorner<3,3>() = vh;
  M.bottomRightCorner<3,3>() = wh;
}

void SE3::adt(Matrix6d &M, const Vector6d &mu) const
{
  Matrix3d mwh;  
  so3.ad(mwh, mu.head<3>());

  Matrix3d mvh;  
  so3.ad(mvh, mu.tail<3>());

  M.topLeftCorner<3,3>() = mwh;
  M.topRightCorner<3,3>() = mvh;
  M.bottomLeftCorner<3,3>() = mvh;
  M.bottomRightCorner<3,3>().setZero();
}



void SE3::adinv(Vector6d& v, const Matrix6d& M) const
{
  Vector3d vw;
  so3.hatinv(vw, M.topLeftCorner<3,3>());

  Vector3d vv;
  so3.hatinv(vv, M.bottomLeftCorner<3,3>());

  v.head<3>() = vw;
  v.tail<3>() = vv;
}


void SE3::exp(Matrix4d &g, const Vector6d &v) const
{
  Matrix3d D;
  so3.dexp(D, v.head<3>());

  Matrix3d R;
  so3.exp(R, v.head<3>());

  g.topLeftCorner<3,3>() = R;
  g.topRightCorner<3,1>() = D*v.tail<3>();
  g.bottomLeftCorner<1,3>().setZero();
  g(3,3) = 1;
}


/*
//f = [atan2(g(2,1), g(1,1)); g(1:2,1:2)'*g(1:2,3)];
void SE3::plog(Vector3d &v, const Matrix3d &g) const
{
  const Vector2d &p = g.block(0,2,2,1);
  v[0] = atan2(g(1,0), g(0,0));
  v.tail(2) = g.topLeftCorner(2,2).transpose()*p;
}
*/


void SE3::log(Vector6d &v, const Matrix4d &g) const
{
  Vector3d w;
  so3.log(w, g.topLeftCorner<3,3>());
  v.head<3>() = w;
  Matrix3d Di;
  so3.dexpinv(Di, w);
  v.tail<3>() = Di*g.topRightCorner<3,1>();
}



void SE3::cay(Matrix4d &g, const Vector6d &v) const
{
  double v0s = v[0]*v[0];
  double v1s = v[1]*v[1];
  double v2s = v[2]*v[2];
  double a = 4 + v0s + v1s + v2s;
  
  double v01 = v[0]*v[1];
  double v02 = v[0]*v[2];
  double v12 = v[1]*v[2];

  g(0,0) = (4+v0s-v1s-v2s)/a;
  g(0,1) = 2*(-2*v[2]+v01)/a;
  g(0,2) = 2*(2*v[1]+v02)/a;
  g(0,3) = (4*v[3]+v[3]*v0s-2*v[4]*v[2]+v[4]*v01+2*v[5]*v[1]+v[5]*v02)/a;
  g(1,0) = 2*(2*v[2]+v01)/a;
  g(1,1) = -(-4+v0s-v1s+v2s)/a;
  g(1,2) = -2*(2*v[0]-v12)/a;
  g(1,3) = (2*v[3]*v[2]+v[3]*v01+4*v[4]+v[4]*v1s-2*v[5]*v[0]+v[5]*v12)/a;
  g(2,0) = 2*(-2*v[1]+v02)/a;
  g(2,1) = 2*(2*v[0]+v12)/a;
  g(2,2) = -(-4+v0s+v1s-v2s)/a;
  g(2,3) = (-2*v[3]*v[1]+v[3]*v02+2*v[4]*v[0]+v[4]*v12+4*v[5]+v[5]*v2s)/a;
  g(3,0) = 0;
  g(3,1) = 0;
  g(3,2) = 0;
  g(3,3) = 1;
}


void SE3::tlnmu(Vector6d& mup, const Vector6d& v, const Vector6d &mu) const
{
  mup.head<3>() = mu.head<3>() + v.head<3>().cross(mu.head<3>())/2;
  mup.tail<3>() = mu.tail<3>() + v.head<3>().cross(mu.tail<3>())/2;
}


void SE3::tln(Matrix6d &M, const Vector6d &v) const
{
  ad(M, v/2);  
  M = Matrix6d::Identity() - M;
}


void SE3::dcayinv(Matrix6d &M, const Vector6d &v) const
{
  Matrix3d whh;
  so3.ad(whh, v.head<3>()/2);

  Matrix3d vhh;  
  so3.ad(vhh, v.tail<3>()/2);

  Matrix3d D;
  so3.dcayinv(D, v.head<3>());

  M.topLeftCorner<3,3>() = D;
  M.topRightCorner<3,3>().setZero();
  M.bottomLeftCorner<3,3>() = -(so3.Id - whh)*vhh;
  M.bottomRightCorner<3,3>() = so3.Id - whh;
}


void SE3::dcay(Matrix6d &M, const Vector6d &v) const
{  
  const double &v0 = v(0);
  const double &v1 = v(1);
  const double &v2 = v(2);
  const double &v3 = v(3);
  const double &v4 = v(4);
  const double &v5 = v(5);

  double n = v0*v1 + v1*v2 + v2*v2 + 4;

  M(0,0) = 4/n;
  M(0,1) -(2*v2)/n;
  M(0,2) = (2*v1)/n;
  M(0,3) = 0;
  M(0,4) = 0;
  M(0,5) = 0;
  M(1,0) = (2*v2)/n;
  M(1,1) = 4/n;
  M(1,2) = -(2*v0)/n;
  M(1,3) = 0;
  M(1,4) = 0;
  M(1,5) = 0;
  M(2,0) = -(2*v1)/n;
  M(2,1) = (2*v0)/n;
  M(2,2) = 4/n;
  M(2,3) = 0;
  M(2,4) = 0;
  M(2,5) = 0;
  M(3,0) = -(v1*v4 + v2*v5)/n;
  M(3,1) = -(2*v5 - v0*v4)/n;
  M(3,2) = (2*v4 + v0*v5)/n;
  M(3,3) = (v0*v0 + 4)/n;
  M(3,4) = -(2*v2 - v0*v1)/n;
  M(3,5) = (2*v1 + v0*v2)/n;
  M(4,0) = (2*v5 + v1*v3)/n;
  M(4,1) = -(v0*v3 + v2*v5)/n;
  M(4,2) = -(2*v3 - v1*v5)/n;
  M(4,3) = (2*v2 + v0*v1)/n;
  M(4,4) = (v1*v1 + 4)/n;
  M(4,5) = -(2*v0 - v1*v2)/n;
  M(5,0) = -(2*v4 - v2*v3)/n;
  M(5,1) = (2*v3 + v2*v4)/n;
  M(5,2) = -(v0*v3 + v1*v4)/n;
  M(5,3) = -(2*v1 - v0*v2)/n;
  M(5,4) = (2*v0 + v1*v2)/n;
  M(5,5) = (v2*v2 + 4)/n;
}

void SE3::q2g(Matrix4d &g, const Vector6d &q) const
{
  Matrix3d R;
  so3.q2g(R, q.head<3>());
  g.topLeftCorner<3,3>() = R;
  g.topRightCorner<3,1>() = q.tail<3>();
  g.bottomLeftCorner<1,3>().setZero();
  g(3,3) = 1;  
}

void SE3::g2q(Vector6d &q, const Matrix4d &g) const
{
  Vector3d o;
  so3.g2q(o, g.topLeftCorner<3,3>());
  q.head<3>() = o;
  q.tail<3>() = g.topRightCorner<3,1>();
}

void SE3::quatxyz2g(Matrix4d &g, const Vector7d &wquatxyz) const
{
	Matrix3d R;
  so3.quat2g(R,wquatxyz.head<4>().normalized());
	g.topLeftCorner<3,3>() = R;
	g.topRightCorner<3,1>() = wquatxyz.tail<3>();
  g.bottomLeftCorner<1,3>().setZero();
  g(3,3) = 1;  
}

void SE3::g2quatxyz(Vector7d &wquatxyz, const Matrix4d &g) const
{
	Vector4d head ;
	so3.g2quat(head,g.topLeftCorner<3,3>());
	wquatxyz.head<4>() = head;
	wquatxyz.tail<3>() = g.topRightCorner<3,1>();
	return;
}


void SE3::rpyxyz2g(Matrix4d &g, const Vector3d &rpy, const Vector3d &xyz) const
{
  Matrix3d R;
  so3.q2g(R, rpy);
  g.topLeftCorner<3,3>() = R;
  g.topRightCorner<3,1>() = xyz;
  g.bottomLeftCorner<1,3>().setZero();
  g(3,3) = 1;  
}

void SE3::g2rpyxyz(Vector3d &rpy, Vector3d &xyz, const Matrix4d &g) const
{
  so3.g2q(rpy, g.topLeftCorner<3,3>());
  xyz = g.topRightCorner<3,1>();
}
