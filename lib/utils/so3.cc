#include "so3.h"
#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace gcop;
using namespace Eigen;
using namespace std;

//vec SO3::e[3] = {vec("1,0,0"), vec("0,1,0"),vec("0,0,1")};
//mat SO3::E[3] = {SO3::ad(e[0]), SO3::ad(e[1]), SO3::ad(e[2])};

SO3::SO3() : Group(),
             tol(1e-16)
{
}

SO3& SO3::Instance()
{
  static SO3 so3;
  return so3;
}



void SO3::inv(Matrix3d &mi, const Matrix3d &m) const
{
  mi = m.transpose();
}


void SO3::Ad(Matrix3d &m, const Matrix3d &g) const
{
  m = g;
}


void SO3::hat(Matrix3d &m, const Vector3d &v) const
{
  m(0,0) = 0; m(0,1) = -v(2); m(0,2) = v(1);
  m(1,0) = v(2); m(1,1) = 0; m(1,2) = -v(0);
  m(2,0) = -v(1); m(2,1) = v(0); m(2,2) = 0;
}


void SO3::hatinv(Vector3d& v, const Matrix3d& m) const
{
  v(0) = m(2,1);
  v(1) = m(0,2);
  v(2) = m(1,0);
}

void SO3::Tg(Matrix3d &m, const Matrix3d &g) const
{
  m = g;
}


void SO3::ad(Matrix3d &m, const Vector3d &v) const
{
  hat(m, v); 
}


void SO3::adinv(Vector3d& v, const Matrix3d& m) const
{
  hatinv(v, m);
}


void SO3::exp(Matrix3d &m, const Vector3d &v) const
{
  double theta = sqrt(v.dot(v));
  if (theta < tol) {
    m.setIdentity();
    return;
  }
  
  Matrix3d vh;
  hat(vh, v);  
  m = Id + sin(theta)/theta*vh + ((1.0-cos(theta))/(theta*theta))*vh*vh;
}


void SO3::log(Vector3d &v, const Matrix3d &m) const
{
  //std::cout << std::setprecision(100) ;
	//cout<<"m :"<<endl<<m<<endl;
	double arg = (m.trace()-1)/2;
	//cout<<" Arg: "<<arg<<endl;
	//assert(abs(arg) <= 1 + tol);
  double phi = 0;
	if(arg >= 1)
		phi = 0;
	else if(arg <= -1)
		phi = M_PI;
	else
		phi = acos(arg);
	//cout<<"phi "<<phi<<endl;
	//cout<<"m trace"<<m.trace()<<endl;
	double sphi = sin(phi);
	//cout<<"sphi "<<sphi<<endl;
  
  if (fabs(sphi) < tol) {
    v.setZero();
    return;    
  }
  
  hatinv(v, (phi/(2.0*sphi))*(m-m.transpose()));
	//cout<<"frac: "<<(phi/(2.0*sphi))<<endl;
}



void SO3::cay(Matrix3d &m, const Vector3d &v) const
{
  Matrix3d vh;
  hat(vh, v);
  m = Id + (1.0/(1+v.dot(v)/4))*(vh + vh*vh/2);
}


void SO3::cayinv(Vector3d &v, const Matrix3d& g) const
{
  // actually this is same as base implementation
  hatinv(v, (Id + g).inverse()*(Id - g));
  v = -2*v;
}


void SO3::dcay(Matrix3d &m, const Vector3d &v) const
{  
  const double& v1 = v[0];
  const double& v2 = v[1];
  const double& v3 = v[2];
  double d = v1*v1 + v2*v2 + v3*v3 + 4;
  
  m(0,0) = 4/d;       m(0,1) = -(2*v3)/d; m(0,2) =  (2*v2)/d;
  m(1,0) = (2*v3)/d;  m(1,1) = 4/d;       m(1,2) = -(2*v1)/d;
  m(2,0) = -(2*v2)/d; m(2,1) = (2*v1)/d;  m(2,2) = 4/d;
}

void SO3::dcayinv(Matrix3d &m, const Vector3d &v) const
{
  const double& v1 = v[0];
  const double& v2 = v[1];
  const double& v3 = v[2];
  m(0,0) = v1*v1/4 + 1;       m(0,1) =  v3/2 + (v1*v2)/4; m(0,2) = (v1*v3)/4 - v2/2;
  m(1,0) = (v1*v2)/4 - v3/2; m(1,1) =  v2*v2/4 + 1;       m(1,2) = v1/2 + (v2*v3)/4;
  m(2,0) = v2/2 + (v1*v3)/4; m(2,1) = (v2*v3)/4 - v1/2; m(2,2) = v3*v3/4 + 1;
//  m = eye(3) - vh/2.0 + itpp::outer_product(v,v)/4.0;
}

void SO3::dexp(Matrix3d &m, const Vector3d &v) const
{  
  double a = sqrt(v.dot(v));
  if (a < tol) {
    m.setIdentity();
    return;
  }

  Matrix3d vh;
  hat(vh, v);

  m = Id + (1-cos(a))*vh/(a*a) + (1-sin(a)/a)*vh*vh/(a*a);
}

void SO3::dexpinv(Matrix3d &m, const Vector3d &v) const
{
  double theta = sqrt(v.dot(v));
  if (theta < tol) {
    m.setIdentity();
    return;
  }
  Matrix3d vh;
  hat(vh, v);
  m = Id - vh/2.0 + ( (1 - theta/(2*tan(theta/2.0))) / (theta*theta))*vh*vh;
}

void SO3::skew(Vector3d& v, const Matrix3d& m) const
{
  hatinv(v, m - m.transpose());
  v = v/2;
}

double SO3::roll(const Matrix3d& m) const
{
  return atan2(m(2,1), m(2,2));
}

double SO3::pitch(const Matrix3d& m) const
{
  return atan2(-m(2,0), sqrt(m(2,1)*m(2,1) + m(2,2)*m(2,2)));
}


double SO3::yaw(const Matrix3d& m) const
{
  return atan2(m(1,0),m(0,0));
}

void SO3::q2g(Matrix3d &m, const Vector3d &rpy)
{
  double ca = cos(rpy(2));
  double sa = sin(rpy(2));
  double cb = cos(rpy(1));
  double sb = sin(rpy(1));
  double cg = cos(rpy(0));
  double sg = sin(rpy(0));
  
  m(0,0) = ca*cb;
  m(0,1) = ca*sb*sg-sa*cg;
  m(0,2) = ca*sb*cg+sa*sg;
  m(1,0) = sa*cb;
  m(1,1) = sa*sb*sg+ca*cg;
  m(1,2) = sa*sb*cg-ca*sg;
  m(2,0) = -sb;
  m(2,1) = cb*sg;
  m(2,2) = cb*cg; 
}
void SO3::g2quat(Vector4d &wxyz, const Matrix3d &m)
{
	//cout<<"M "<<endl<<m<<endl;
	double trace = m(0,0) + m(1,1) + m(2,2);
	double temp[4];

	if (trace > double(0.0)) 
	{
		double s = sqrt(trace + double(1.0));
		temp[3]=(s * double(0.5));
		s = double(0.5) / s;

		temp[0]=((m(2,1) - m(1,2)) * s);
		temp[1]=((m(0,2) - m(2,0)) * s);
		temp[2]=((m(1,0) - m(0,1)) * s);
	} 
	else 
	{
		int i = m(0,0) < m(1,1) ? 
			(m(1,1) < m(2,2) ? 2 : 1) :
			(m(0,0) < m(2,2) ? 2 : 0); 
		int j = (i + 1) % 3;  
		int k = (i + 2) % 3;

		double s = sqrt(m(i,i) - m(j,j) - m(k,k) + double(1.0));
		//cout<<"s "<<s<<endl;
		temp[i] = s * double(0.5);
		s = double(0.5) / s;

		temp[3] = (m(k,j) - m(j,k)) * s;
		temp[j] = (m(j,i) + m(i,j)) * s;
		temp[k] = (m(k,i) + m(i,k)) * s;
	}
	wxyz<<temp[3],temp[0],temp[1],temp[2];
	//cout<<"wxyz "<<wxyz<<endl;
	//q.setValue(temp[0],temp[1],temp[2],temp[3]);



/*
	if (trace > 0.0) 
	{
		double s = sqrt(trace + 1.0);
		wxyz(0)=(s * double(0.5));
		s = double(0.5) / s;

		wxyz(1)=((m(2,1) - m(1,2)) * s);
		wxyz(2)=((m(0,2) - m(2,0)) * s);
		wxyz(3)=((m(1,0) - m(0,1)) * s);
	} 
	else 
	{
		int i = m(0,0) < m(1,1) ? 
			(m(1,1) < m(2,2) ? 2 : 1) :
			(m(0,0) < m(2,2) ? 2 : 0); 
		int j = (i + 1) % 3;  
		int k = (i + 2) % 3;

		double s = sqrt(m(i,i) - m(j,j) - m(k,k) + double(1.0));
		wxyz(i+1) = s * double(0.5);
		s = double(0.5) / s;

		wxyz(0) = (m(k,j) - m(j,k)) * s;
		wxyz(j) = (m(j,i) + m(i,j)) * s;
		wxyz(k) = (m(k,i) + m(i,k)) * s;
	}
	//q.setValue(temp[0],temp[1],temp[2],temp[3]);
	*/
}
void SO3::quat2g(Matrix3d &m, const Vector4d &wxyz)
{
	assert(abs(wxyz.norm() - 1) < 1e-5);
	double xy = wxyz(1)*wxyz(2), yz = wxyz(3)*wxyz(2),wy =wxyz(2)*wxyz(0),wx = wxyz(1)*wxyz(0),wz = wxyz(0)*wxyz(3), xz = wxyz(1)*wxyz(3), x2 = wxyz(1)*wxyz(1), y2 = wxyz(2)*wxyz(2), z2 = wxyz(3)*wxyz(3), w2 = wxyz(0)*wxyz(0);
	m(0,0) = 1 - 2*(y2 + z2);
  m(0,1) = 2*(xy - wz);
  m(0,2) = 2*(xz + wy);
  m(1,0) = 2*(xy + wz);
  m(1,1) = 1 - 2*(x2 + z2);
  m(1,2) = 2*(yz - wx);
  m(2,0) = 2*(xz - wy);
  m(2,1) = 2*(yz + wx);
  m(2,2) = 1 - 2*(x2 + y2); 
}

void SO3::g2q(Vector3d &rpy, const Matrix3d &m)
{
  rpy[0] = SO3::roll(m);
  rpy[1] = SO3::pitch(m);
  rpy[2] = SO3::yaw(m);
}
