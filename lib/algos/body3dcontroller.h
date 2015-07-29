#ifndef GCOP_BODY3DCONTROLLER_H
#define GCOP_BODY3DCONTROLLER_H

#include "controller.h"
#include "body3d.h"
#include "so3.h"

namespace gcop {

  using namespace std;
  using namespace Eigen;

  /**
   * Rigid body PD controller using errors on \f$ SO(3) \times R^3\f$
   * 
   * Author: Marin Kobilarov, 2007
   */
  template <int c = 6> 
    class Body3dController : public Controller<Body3dState, Matrix<double, c, 1> > {
  public:
  typedef Matrix<double, c, 1> Vectorcd;
  
  /**
   * BODY3D PD controller
   *
   * @param sys multi-body system
   * @param xd desired state (optional, set to origin by default)
   * @param ad desired acceleration (optional, set to zero by default)
   */
  Body3dController(const Body3d<c> &sys,
                   Body3dState *xd = 0, 
                   Vector6d *ad = 0);

  virtual bool Set(Vectorcd &u, double t, const Body3dState &x);
  
  virtual ~Body3dController();
  
  const Body3d<c> &sys; ///< system
  
  Body3dState *xd;   ///< desired state (origin by default)
  Vector6d *ad;      ///< desired acceleration (zero by default)
  
  Vector6d Kp;    ///< proportional terms  (ones by default)
  Vector6d Kd;    ///< derivative  terms  (twos by default)
  
  };
  
  template <int c>   
    Body3dController<c>::Body3dController(const Body3d<c> &sys,
                                          Body3dState *xd, 
                                          Vector6d *ad) :
    sys(sys), xd(xd), ad(ad)
  {
    Kp << 1, 1, 1, 1, 1, 1;
    Kd << 2, 2, 2, 2, 2, 2;
  }
  
  template <int c>   
    Body3dController<c>::~Body3dController()
    {      
    }
  
  template <int c>   
    bool Body3dController<c>::Body3dController::Set(Vectorcd &u, double t, const Body3dState &x)
    {
      const Matrix3d& R = x.R;  // rotation

      Vector3d eR;  // error in rotation
      Vector3d ew;  // error in ang vel (body frame)
      Vector3d ex;  // error in position
      Vector3d ev;  // error in lin vel (spatial frame)

      if (xd) {
        SO3::Instance().log(eR, xd->R.transpose()*R);
        ew = x.w - x.R.transpose()*xd->R*xd->w;
        ex = x.p - xd->p;
        ev = x.v - xd->v;
      } else {
        SO3::Instance().log(eR, R);
        ew = x.w;
        ex = x.p;
        ev = x.v;
      }

      u.head(3) = -Kp.head<3>().cwiseProduct(eR) - Kd.head<3>().cwiseProduct(ew);
      u.tail(3) = R.transpose()*(-Kp.tail<3>().cwiseProduct(ex) - Kd.tail<3>().cwiseProduct(ev) - sys.fp);
      return true;
    }
};

#endif
