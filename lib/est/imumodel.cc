#include "imumodel.h"
#include "quat.h"
#include <assert.h>

using namespace gcop;
using namespace std;
using namespace Eigen;


ImuModel::ImuModel(bool acc, bool gmr) :
  Model((acc ? 10 : 7), (acc ? 6 : 3), 3, 1),
  acc(acc), gmr(gmr)
{
  quat = true;

  
// various noise configurations

#define IMU_HIP

#ifdef IMU_HIP
  
  // output noise : 0.9 deg/s= 0.0157 rad/s
  wn <<.0157 ,.0157 ,.0157;

  // initial sensitivity  : 0.05 deg/s= 0.0157 rad/s
  //  wn0 = ".0157 .0157 .0157";

  // bias stability: 
  //  wbn = ".0005 .0005 .0005";
  wbn << .05, .05, .05;

  // initial bias error : 3 deg/s = 0.0524
  //  wbn0 = ".05 .05 .05";
  wbn0 << .05, .05, .05;
  
  // output noise: 9 mg = 0.082 m/s^2
  vn << .08, .08, .08;

  // acc bias stability: 0.2mg = .002 m/s^2
  //  abn = ".0002 .0002 .0002";
  abn << .02, .02, .02;

  // initial bias error: 50mg = .49 m/s^2
  abn0 << .5, .5, .5;

  // mag bias stability
  //  mbn = ".0084 .0084 .0084";

  // initial bias error 4 mg = 4/475 = .0084
  //  mbn0 = ".0084 .0084 .0084";

  // output noise 1.25 mg = 1.25/475
  mon << .0026, .0026, .0026;

#else
  wn = "0.5 0.5 0.5";
  wbn ="0.05 0.05 0.05";
  vn = "1 1 1";
  abn = ".01 .01 .01";
  mon = ".1 .1 .1";

  wbn0 = ".5 .5 .5";
  abn0 = ".5 .5 .5";

#endif

  // gravity
  g << 0, 0, -9.81;

  // magnetic north reference
  // in nT  (1e9 nT  = 1e4 G )
  //  m0 = "23984.0 5411.0 -40647.0";   // nT 
  m0 << 0.23984, 0.05411, 0.40647;   // G

  // total field 47,505.6  nT (0.475 G)
  m0.normalize();

}

ImuModel::~ImuModel()
{
}


bool ImuModel::f(VectorXd &xn,
                 const VectorXd &x,
                 const VectorXd &u)
{
  assert(dt > 0);

  Quat q(x.data());
  Vector3d wb(x.data() + 4);  // gcop gyro bias
  Vector3d wo(u.data());      // observed gyro
  Vector3d w = wo - wb;   // corrected ang. vel.
  
  Quat dq;
  Vector3d dtw = this->dt*w;
  dq.FromExp(dtw.data());
  Quat qn = q*dq;
  qn.Normalize();

  qn.Q(xn.data());
  memcpy(xn.data() + 4, x.data() + 4, 3*sizeof(double));   // bias is the same

  if (acc) {
    memcpy(xn.data() + 7, x.data() + 7, 3*sizeof(double)); // acc is the same
  }
  return true;
}


bool ImuModel::h(VectorXd &z, 
                 const VectorXd &x)
{
  Quat qi(x.data());  
  qi.Invert();

  qi.Rotate2(z.data(), m0.data());   // mo = q2R(q)'*m0

  if (acc) {
    Vector3d ab(x.data() + 7);  // acc bias
    Vector3d a;          // corrected acc  
    qi.Rotate2(a.data(), g.data());  
    // ao = R'*a0 + ab
    z[3] = a[0] + ab[0];
    z[4] = a[1] + ab[1];
    z[5] = a[2] + ab[2];
  }
  return true;
}

void ImuModel::I(double *y,
                 const double *z,
                 const double *x)
{
  Quat q(x);  
  q.Invert();

  double m[3];
  q.Rotate2(m, m0.data());   // mo = q2R(q)'*m0

  y[0] = z[0] - m[0];
  y[1] = z[1] - m[1];
  y[2] = z[2] - m[2];

  if (acc) {
    Vector3d ab(x + 7);  // acc bias
    Vector3d a(3);          // corrected acc
    q.Rotate2(a.data(), g.data());
    // ao = R'*a0 + ab
    y[3] = z[3] - (a[0] + ab[0]);
    y[4] = z[4] - (a[1] + ab[1]);
    y[5] = z[5] - (a[2] + ab[2]);
  }
}

bool ImuModel::IsValid(const Vector3d &w,
                       const Vector3d &a,
                       const Vector3d &m) const
{
  for (int i = 0; i < 3; ++i) {
    if (w.size() && fabs(w[i]) > 300*M_PI/180.0)
      return false;
    if (a.size() && fabs(a[i]) > 18*9.81)
      return false;
    if (m.size() && fabs(m[i]) > 1)
      return false;    
  }
  return true;
}
