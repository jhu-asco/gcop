#include "imuukf.h"
#include <assert.h>
#include <quat.h>

using namespace gcop;
using namespace std;
using namespace Eigen;

ImuUKF::ImuUKF(ImuModel &im) : UKF(im), z(6)
{
  t = -1;
  Reset();
}

ImuUKF::~ImuUKF()
{
}

void ImuUKF::Reset()
{
  ImuModel& im = (ImuModel&)this->model;

  x.setZero();
  x[0] = 1;
  
  if (im.acc) {
    P.diagonal() << .001*im.wn.cwiseProduct(im.wn), im.wbn0.cwiseProduct(im.wbn0), im.abn0.cwiseProduct(im.abn0);  
  } else {
    P.diagonal() << .001*im.wn.cwiseProduct(im.wn), im.wbn0.cwiseProduct(im.wbn0);
  }
  updated = false;
}


bool ImuUKF::Process(double t,
                     const Vector3d &w,
                     const Vector3d &a,
                     const Vector3d &m)  
{
  ImuModel& im = (ImuModel&)this->model;
  
  if (!im.IsValid(w, a, m)) {
    cout << "[W] ImuUkf::Process: (t=" << t << ")-- measurement not valid!" << endl;
    return false;
  }
  
  double dt = t - this->t;
  this->t = t;
  
  if (!updated) {
    double *x_ = this->x.data();
    
    if (!im.gmr) {
      im.m0[0] = m[0];
      im.m0[1] = m[1];
      im.m0[2] = m[2];
    }
    
    memcpy(x_ + 4, w.data(), 3*sizeof(double));  // gyro
    
    updated = true;
    return true;
  }
  
  assert(dt > 0);

  // update model (i.e. dt, Q)
  //  im.Update(dt, false);

  im.dt = dt;

  if (im.acc) {
    im.Q.diagonal() << dt*im.wn.cwiseProduct(im.wn) + dt*dt*dt/3*im.wbn.cwiseProduct(im.wbn), 
      dt*im.wbn.cwiseProduct(im.wbn), dt*im.abn.cwiseProduct(im.abn);
    
  } else {
    im.Q.diagonal() << dt*im.wn.cwiseProduct(im.wn), 
      dt*im.wbn.cwiseProduct(im.wbn);
  }
  // cross-terms (not very critical)
  im.Q.block(0,3,3,3).diagonal() << -dt*dt/2*im.wbn.cwiseProduct(im.wbn);
  im.Q.block(3,0,3,3).diagonal() << -dt*dt/2*im.wbn.cwiseProduct(im.wbn);
  
  
  // predict
  VectorXd u(w);
  Predict(&u);

  if (im.acc) {
    im.R.diagonal() << im.mon.cwiseProduct(im.mon), im.vn.cwiseProduct(im.vn);
  } else {
    im.R.diagonal() << im.mon.cwiseProduct(im.mon);
  }
  
  memcpy(z.data(), m.data(), 3*sizeof(double));
  if (im.acc)
    memcpy(z.data() + 3, a.data(), 3*sizeof(double));
  
  // update
  return Update(z);
}
