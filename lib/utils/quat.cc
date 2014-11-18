#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include "quat.h"
#include <limits>

//#define EPS (numeric_limits<double>::epsilon())
#define EPS 1e-16

using namespace std;
using namespace gcop;

Quat::Quat()  
{ 
  Identity(); 
  tol = EPS;
}

Quat::Quat(const Quat& q) : qw(q.qw), qx(q.qx), qy(q.qy), qz(q.qz), tol(EPS)
{
}

Quat::Quat(const double *q) : qw(q[0]),
                              qx(q[1]),
                              qy(q[2]), 
                              qz(q[3]),
                              tol(EPS)
{
  //  Normalize();
}

Quat::Quat(double w, double x, double y, double z) : qw(w),
                                                     qx(x), 
						     qy(y), 
						     qz(z),
                                                     tol(EPS)
{
}

Quat::~Quat()
{
}

void Quat::Identity() 
{
  qw = 1;
  qx = qy = qz = 0;
}

double* Quat::Rotate(double c[3], const double a[3]) const
{
  if (IsIdentity()) {
    if (a != c)
      memcpy(c, a, 3*sizeof(double));
    return c;
  }
 
  double n = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  if (n <= tol) {
    cerr << "[W]: Quat::Rotate: vector norm close to zero" << endl;
    memcpy(c, a, 3*sizeof(double));
    return c;    
  }

  Quat q(0, a[0]/n, a[1]/n, a[2]/n);
  Quat i; Invert(i);
  q = (*this)*q*i;
  q.ToExp(c);
  
  c[0] *= n;
  c[1] *= n;
  c[2] *= n;

  return c;
}


double* Quat::Rotate(double a[3]) const
{
  return Rotate(a, a);
}


double* Quat::Rotate(double c[3], const double a[3], const double u[3], double v)
{
  Quat q;
  q.FromAxis(u, v);
  return q.Rotate(c, a);
}

double* Quat::Rotate2(double c[3], const double a[3]) const
{
  c[0] = (qw*qw + qx*qx - qy*qy -qz*qz)*a[0] + (2*(qx*qy - qw*qz))*a[1] + (2*(qx*qz + qw*qy))*a[2];
  c[1] = (2*(qx*qy + qw*qz))*a[0] + (qw*qw - qx*qx + qy*qy - qz*qz)*a[1] + (2*(qy*qz - qw*qx))*a[2];
  c[2] = (2*(qx*qz - qw*qy))*a[0] + (2*(qy*qz + qw*qx))*a[1] + (qw*qw - qx*qx - qy*qy + qz*qz)*a[2];
  return c;
}


double* Quat::Rotate2(double c[3]) const
{
  double a[3];
  a[0] = c[0];
  a[1] = c[1];
  a[2] = c[2];
  return Rotate2(c, a);
}

void Quat::Transform(const double e[3]) 
{
  Quat q;
  q.FromExp(e);
  *this = q*(*this);
}


void Quat::Transform(const double u[3], const double v) 
{
  Quat q;
  q.FromAxis(u, v);
  *this = q*(*this);
}
 
void Quat::Invert() 
{
  qx = -qx;
  qy = -qy;
  qz = -qz;
}

void Quat::Invert(Quat &q) const
{
  q.qw = qw;
  q.qx = -qx;
  q.qy = -qy;
  q.qz = -qz;
}

void Quat::Normalize()
{
  double norm = sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
  qw /= norm;
  qx /= norm;
  qy /= norm;
  qz /= norm;
}


Quat& Quat::operator*=(const Quat &q)
{
  if (q.IsIdentity())
    return *this;
  double t0 = (qz-qy)*(q.qy-q.qz);
  double t1 = (qw+qx)*(q.qw+q.qx);
  double t2 = (qw-qx)*(q.qy+q.qz);
  double t3 = (qz+qy)*(q.qw-q.qx);
  double t4 = (qz-qx)*(q.qx-q.qy);
  double t5 = (qz+qx)*(q.qx+q.qy);
  double t6 = (qw+qy)*(q.qw-q.qz);
  double t7 = (qw-qy)*(q.qw+q.qz);
  double t8 = t5+t6+t7;
  double t9 = (t4+t8)/2.0;
  qw = t0+t9-t5;
  qx = t1+t9-t8;
  qy = t2+t9-t7;
  qz = t3+t9-t6;
  
  //  Normalize();
  return *this;
}


void Quat::FromExp(const double e[3]) 
{
  double n = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
  qw = cos(n/2);
  double sn = sin(n/2)/n;
  qx = e[0]*sn;
  qy = e[1]*sn;
  qz = e[2]*sn;
  //  Normalize();
}

void Quat::ToExp(double e[3]) const
{
  if (IsIdentity()) {
    e[0] = e[1] = e[2] = 0;    
    return;
  }
  
  double ac = acos(qw);
  if (ac < tol) {
    e[0] = 0;
    e[1] = 0;
    e[2] = 0;
    return;
  }

  double sn = sin(ac)/(2*ac);
  e[0] = qx/sn;
  e[1] = qy/sn;
  e[2] = qz/sn;
}



void Quat::FromAxis(const double u[3], const double v) 
{
  qw = cos(v/2);
  double sinv = sin(v/2);
  qx = u[0]*sinv;
  qy = u[1]*sinv;
  qz = u[2]*sinv;
  Normalize();
}


void Quat::ToAxis(double u[3], double &v) const
{
  if (IsIdentity()) {
    v = 0;
    if (u) {
      u[0] = 0; u[1] = 0; u[2] = 0;
    }
    return;
  }
  
  double ac = acos(qw);
  v = 2*ac;

  double sinv = sin(ac);
  u[0] = qx/sinv;
  u[1] = qy/sinv;
  u[2] = qz/sinv;
}



 
void Quat::ToSE3(double m[16]) const
{
  m[0] =  (qw*qw + qx*qx - qy*qy -qz*qz);
  m[1] =  (2*(qx*qy + qw*qz));
  m[2] =  (2*(qx*qz - qw*qy));
 
  m[4] =  (2*(qx*qy - qw*qz));
  m[5] =  (qw*qw - qx*qx + qy*qy - qz*qz);
  m[6] =  (2*(qy*qz + qw*qx));
   
  m[8] =  (2*(qx*qz + qw*qy));
  m[9] =  (2*(qy*qz - qw*qx));
  m[10]=  (qw*qw - qx*qx - qy*qy + qz*qz);
   
  m[15]=  (qw*qw + qx*qx + qy*qy + qz*qz);
   
  m[3] = m[7] = m[11] = m[12] = m[13] = m[14] = 0;
}
 
static int perm[3] = {1,2,0};
#define mat(a,b) (m[4*a+b])
 
void Quat::FromSE3(const double m[16]) 
{
  double T = mat(0,0) + mat(1,1) + mat(2,2);
  if (T > 0) { // w is the largest element in the quat
    double iw = sqrt(T+1.0);
    qw = iw * 0.5; 
    iw = 0.5/iw;
    qx = (mat(1,2) - mat(2,1)) * iw;
    qy = (mat(2,0) - mat(0,2)) * iw;
    qz = (mat(0,1) - mat(1,0)) * iw;
  } else {     // Find the largest diagonal element
    int i,j,k;
    double &qi = qx, &qj = qy, &qk = qz;
    i=0;
    if (mat(1,1) > mat(0,0)) {i=1; qi = qy; qj = qz; qk = qx; }
    if (mat(2,2) > mat(i,i)) {i=2; qi = qz; qj = qx; qk = qy; }
    j=perm[i];
    k=perm[j];
  
    double iqi = sqrt( (mat(i,i) - (mat(j,j) + mat(k,k))) + 1.0);
    qi = iqi * 0.5;
    iqi = 0.5/iqi;
    
    qw = (mat(j,k) - mat(k,j))*iqi;
    qj = (mat(i,j) + mat(j,i))*iqi;
    qk = (mat(i,k) + mat(k,i))*iqi;
  }
  Normalize();
}



void Quat::FromRpy(const double* rpy)
{
  double hr = rpy[0]/2;
  double hp = rpy[1]/2;
  double hy = rpy[2]/2;
  double chr = cos(hr);
  double shr = sin(hr);
  double chp = cos(hp);
  double shp = sin(hp);
  double chy = cos(hy);
  double shy = sin(hy);
  
  qw = chr*chp*chy + shr*shp*shy;
  qx = shr*chp*chy - chr*shp*shy;
  qy = chr*shp*chy + shr*chp*shy;
  qz = chr*chp*shy - shr*shp*chy;
  Normalize();
}

void Quat::ToRpy(double* rpy) const
{
  if (IsIdentity()) {
    rpy[0] = rpy[1] = rpy[2] = 0;
    return;
  }
  rpy[0] = atan2(2*(qw*qx + qy*qz), 1 - 2*(qx*qx + qy*qy));
  rpy[1] = asin(2*(qw*qy - qx*qz));
  rpy[2] = atan2(2*(qx*qy + qw*qz), 1 - 2*(qy*qy +qz*qz));
}


double* Quat::Q(double *q) const
{ 
  q[0] = qw; q[1] = qx; q[2] = qy; q[3]= qz; 
  return q;
}

bool Quat::IsIdentity() const
{
  //  return qx*qx+qy*qy+qz*qz <= tol && fabs(qw-1) <= tol;
  //return fabs(qw-1) <= tol;
  return (fabs(qx) < tol && fabs(qy) < tol && fabs(qz) < tol) || fabs(fabs(qw)-1) < tol;
}

void Quat::SetTol(double tol)
{
  this->tol = tol;
}

namespace gcop {

  const Quat operator*(const Quat& q0, const Quat& q1)
  {
    Quat q2 = q0;
    return q2 *= q1;
  }

  std::ostream& operator<<(std::ostream &os, const Quat &q)
  {
    os << q.qw << " " << q.qx << " " << q.qy << " " << q.qz;
    return os;
  }
  
  std::istream& operator>>(std::istream &is, Quat &q)
  {
    is >> q.qw;
    is >> q.qx;
    is >> q.qy;
    is >> q.qz;
    return is;
  }
}
