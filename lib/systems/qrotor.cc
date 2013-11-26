#include "qrotor.h"

using namespace gcop;
using namespace Eigen;

Qrotor::Qrotor() : 
  Body3d(Vector3d(.5, .5, .1), .5, Vector3d(2.32e-3, 2.32e-3, 4.41e-3)),
  l(0.171), r(0.099), 
  kt(6.107e-8), km(2.7207e-9), mm(20) {
  
  Bu.setZero();
  Bu(0,0) = 1;
  Bu(1,1) = 1;
  Bu(2,2) = 1;
  Bu(5,3) = 1;
  
  fp << 0, 0, -9.81*m;

  //  double ulb[4] = {1200, 1200, 1200, 1200};
  //  double uub[4] = {7800, 7800, 7800, 7800};
}


/*


void Qrotor::Control::Bu(double *f, const State &s)
{
  const Qrotor &qrotor = (const Qrotor&)s.sys;
  
  // forces produced by each rotor
  double rfs[4];
  rfs[0] = qrotor.kt*s.u[0]*s.u[0];
  rfs[1] = qrotor.kt*s.u[1]*s.u[1];
  rfs[2] = qrotor.kt*s.u[2]*s.u[2];
  rfs[3] = qrotor.kt*s.u[3]*s.u[3];

  f[0] = qrotor.l*(rfs[3] - rfs[1]);
  f[1] = qrotor.l*(rfs[2] - rfs[0]);
  //  f[2] = qrotor.km*(s.u[0]*s.u[0] - s.u[1]*s.u[1] + s.u[2]*s.u[2] - s.u[3]*s.u[3]);
  f[2] = qrotor.km/qrotor.kt*(rfs[0] - rfs[1] + rfs[2] - rfs[3]);
  f[3] = 0;
  f[4] = 0;
  f[5] = rfs[0] + rfs[1] + rfs[2] + rfs[3]; 
}
*/
