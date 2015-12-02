#include "se2.h"
#include <assert.h>

using namespace gcop;
using namespace Eigen;
using namespace std;

//vec SE2::e[3] = {vec("1,0,0"), vec("0,1,0"),vec("0,0,1")};
//mat SE2::E[3] = {SE2::ad(e[0]), SE2::ad(e[1]), SE2::ad(e[2])};

SE2::SE2() : Group(),
             tol(1e-16)
{
}

SE2& SE2::Instance()
{
  static SE2 se2;
  return se2;
}



void SE2::inv(Matrix3d &mi, const Matrix3d &m) const
{
  const double& ct = m(0,0);
  const double& st = m(1,0);
  const double& x = m(0,2);
  const double& y = m(1,2);
  
  mi(0,0) = ct;   mi(0,1) = st;   mi(0,2) = -x*ct-y*st;
  mi(1,0) = -st;  mi(1,1) = ct;   mi(1,2) = x*st - y*ct;
  mi(2,0) = 0;    mi(2,1) = 0;    mi(2,2) = 1;
}


void SE2::q2g(Matrix3d &m, const Vector3d &q) const
{
  double ct = cos(q[0]);
  double st = sin(q[0]);

  m(0,0) = ct;   m(0,1) = -st;  m(0,2) = q[1];
  m(1,0) = st;   m(1,1) = ct;   m(1,2) = q[2];
  m(2,0) = 0;    m(2,1) = 0;    m(2,2) = 1;
}

void SE2::g2q(Vector3d &q, const Matrix3d &m) const
{
  q[0] = atan2(m(1,0), m(0,0));
  q[1] = m(0,2);
  q[2] = m(1,2);
}

void SE2::Ad(Matrix3d &m, const Matrix3d &g) const
{
  const double &c = g(0,0);
  const double &s = g(1,0);
  const double &x = g(0,2);
  const double &y = g(1,2);
  
  m(0,0) = 1;
  m(0,1) = m(0,2) = 0;
  m(1,0) = y;
  m(1,1) = c;
  m(1,2) = -s;
  m(2,0) = -x;
  m(2,1) = s;
  m(2,2) = c;
}


void SE2::hat(Matrix3d &m, const Vector3d &v) const
{
  m(0,0) = 0; m(0,1) = -v(0); m(0,2) = v(1);
  m(1,0) = v(0); m(1,1) = 0; m(1,2) = v(2);
  m(2,0) = 0; m(2,1) = 0; m(2,2) = 0;
}


void SE2::hatinv(Vector3d& v, const Matrix3d& m) const
{
  v(0) = m(1,0);
  v(1) = m(0,2);
  v(2) = m(1,2);
}

void SE2::Tg(Matrix3d &m, const Matrix3d &g) const
{
  m.bottomRightCorner(2,2) = g.topLeftCorner(2,2);
  m(0,0) = 1; m(0,1) = 0; m(0,2) = 0;
  m(1,0) = 0; 
  m(2,0) = 0;   
}


void SE2::ad(Matrix3d &m, const Vector3d &v) const
{
  m(0,0) = 0;     m(0,1) = 0;    m(0,2) = 0;
  m(1,0) = v(2);  m(1,1) = 0;    m(1,2) = -v(0);
  m(2,0) = -v(1); m(2,1) = v(0); m(2,2) = 0;
}


void SE2::adinv(Vector3d& v, const Matrix3d& m) const
{
  v(0) = m(2,1);
  v(1) = -m(2,0);
  v(2) = m(1,0);
}


void SE2::exp(Matrix3d &m, const Vector3d &v) const
{
  const double &w = v[0];

  if (fabs(w) < tol) {
    m(0,0) = 1; m(0,1) = 0; m(0,2) = v(1);
    m(1,0) = 0; m(1,1) = 1; m(1,2) = v(2);
    m(2,0) = 0; m(2,1) = 0; m(2,2) = 1;
    return;
  }
  
  double c = cos(w);
  double s = sin(w);
  
  double ax = v[2]/w;
  double ay = -v[1]/w;

  m(0,0) = c; m(0,1) = -s; m(0,2) = (c - 1)*ax - s*ay;
  m(1,0) = s; m(1,1) = c;  m(1,2) = s*ax + (c - 1)*ay;
  m(2,0) = 0; m(2,1) = 0;  m(2,2) = 1;
}

//f = [atan2(g(2,1), g(1,1)); g(1:2,1:2)'*g(1:2,3)];
void SE2::plog(Vector3d &v, const Matrix3d &g) const
{
  const Vector2d &p = g.block(0,2,2,1);
  v[0] = atan2(g(1,0), g(0,0));
  v.tail(2) = g.topLeftCorner(2,2).transpose()*p;
}



void SE2::log(Vector3d &v, const Matrix3d &m) const
{
  double th = atan2(m(1,0),m(0,0));
  v(0) = th;
  const double& x = m(0,2);
  const double& y = m(1,2);
  if (fabs(th) < tol) {
    v(1) = x;
    v(2) = y;
    return;
  }
  double th2 = th/2;
  double a = th2/tan(th/2);
  v(1) = a*x + th2*y;
  v(2) = -th2*x + a*y;
}



void SE2::cay(Matrix3d &m, const Vector3d &v) const
{
  const double& w = v(0);
  const double& vx = v(1);
  const double& vy = v(2);
  double t = 1/(4 + w*w);
  double ct = t*(4 - w*w);
  double st = t*4*w;
  m(0,0) = ct;
  m(1,0) = st;
  m(2,0) = 0;
  m(0,1) = -st;
  m(1,1) = ct;
  m(2,1) = 0;
  m(0,2) = -t*2*(w*vy - 2*vx);
  m(1,2) = t*2*(2*vy + vx*w);
  m(2,2) = 1;
}

void SE2::cayinv(Vector3d &v, const Matrix3d& g) const
{
  const double &c = g(0,0);
  const double &s = g(1,0);
  const double &x = g(0,2);
  const double &y = g(1,2);

  if (fabs(c + 1) < tol) {
    v[0] = 0;
    v[1] = x;
    v[2] = y;
  } else {
    double a = s/(1 + c);
    v[0] = 2*a;
    v[1] = x + a*y; 
    v[2] = y - a*x;
  }
}


void SE2::dcayinv(Matrix3d &m, const Vector3d &v) const
{
  //  mat vh = ad(v);
  //  m = (1+v*v/4)*eye(3) - vh/2 + vh*vh/4;
  const double& w = v[0];
  const double& vx = v[1];
  const double& vy = v[2];

  m(0,0) = 1+w*w/4;
  m(1,1) = 1; m(2,2) = 1;
  m(0,1) = 0; m(0,2) = 0;
  m(1,0) = -vy/2 + w*vx/4;
  m(1,2) = w/2;
  m(2,0) = vx/2 + w*vy/4;
  m(2,1) = -w/2;
}

void SE2::dcay(Matrix3d &m, const Vector3d &v) const
{  
  const double &w = v[0];
  double d = 4 + w*w;
  m(0,0) = 4/d;       m(0,1) = 0;      m(0,2) = 0;
  m(1,0) = 2*v[2]/d;  m(1,1) = 4/d;    m(1,2) = -2*w/d;
  m(2,0) = -2*v[1]/d; m(2,1) = 2*w/d;  m(2,2) = 4/d;
}
