#ifndef QROTORIDMODEL_H
#define QROTORIDMODEL_H

#include "qrotoridmanifold.h"
#include "body3d.h"
#include "so3.h"

namespace gcop {

  using namespace std;
  using namespace Eigen;

  class QRotorIDModel: public System<QRotorIDState, 15, 4,10> {
    protected:
      typedef Matrix<double, 4, 1> Vectorcd;
      typedef Matrix<double, 10, 1> Vectormd;
      typedef Matrix<double, 15, 15> Matrixnd;
      typedef Matrix<double, 15, 4> Matrixncd;
      typedef Matrix<double, 15, 10> Matrixnmd;
  protected:
      Vector3d e3;///< Axis direction of quadrotor
      Vector3d g;///< Gravity Vector
      SO3 &so3;///< Instance of so3
    public:
      double kt;///< Thrust gain for quadrotor
      Vector3d kp;///< Proportional gains for rpy control
      Vector3d kd;///< Derivative gains for rpy control
      Vector3d a0;///< Offsets in global acceleration due to wind etc (Add offsets in torque later)
      //, kp(5,5,5)
      //, kd(7,7,7)
    public:
      QRotorIDModel():System<QRotorIDState,15,4,10>(QRotorIDManifold::Instance())
                      , so3(SO3::Instance())
                      , kt(0.16)
                      , kp(6.76,7.14,7.27)
                      , kd(5.65,5.62,5.62)
                      , e3(0,0,1)
                      , a0(0,0,0)
                      , g(0,0,-9.81)
      {
      }
      double Step(QRotorIDState &xb, double t, const QRotorIDState &xa, const Vectorcd &u, double h, const Vectormd *p, Matrixnd *A, Matrixncd *B, Matrixnmd *C)
      {
          Vector3d rpy;
          so3.g2q(rpy,xa.R);
          Vector3d erpy = rpy - xa.u;
          Vector3d eomega = xa.w - u.tail<3>();
          if(p == 0)//No Parameters provided
          {
              //Translational Part
              xb.v = xa.v + h*((kt*u(0))*(xa.R * e3) + g + a0);
              //Rotational Part
              Vector3d omegadot = -kp.cwiseProduct(erpy) - kd.cwiseProduct(eomega);
              xb.w = xa.w + omegadot*h;
          }
          else
          {
            //Params p is 10x1 vector with kt, kp, kd, a0
            //Translational Part
            xb.v = xa.v + h*(((*p)(0)*u(0))*(xa.R * e3) + g + p->tail<3>());
            //Rotational Part
            Vector3d omegadot = -p->segment<3>(1).cwiseProduct(erpy) - p->segment<3>(4).cwiseProduct(eomega);
            xb.w = xa.w + omegadot*h;
          }
          xb.p = xa.p + 0.5*(xb.v + xa.v)*h;
          Vector3d w_avg = 0.5*(xb.w + xa.w);
          Matrix3d dR;
          so3.exp(dR,w_avg*h);
          xb.R = xa.R*dR;
          xb.u = xa.u + u.tail<3>()*h;
          //Specify A,B,C also TODO
      }
  };
}

#endif // QROTORIDMODEL_H
