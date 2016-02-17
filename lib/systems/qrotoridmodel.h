#ifndef QROTORIDMODEL_H
#define QROTORIDMODEL_H

#include "qrotoridmanifold.h"
#include "body3d.h"
#include "so3.h"

namespace gcop {

  using namespace std;
  using namespace Eigen;

  class QRotorIDModel: public System<QRotorIDState, 15, 4,13> {
    protected:
      typedef Matrix<double, 4, 1> Vectorcd;
      typedef Matrix<double, 13, 1> Vectormd;
      typedef Matrix<double, 15, 15> Matrixnd;
      typedef Matrix<double, 15, 4> Matrixncd;
      typedef Matrix<double, 15, 13> Matrixnmd;
  protected:
      Vector3d e3;///< Axis direction of quadrotor
      Vector3d g;///< Gravity Vector
      SO3 &so3;///< Instance of so3
    public:
      double kt;///< Thrust gain for quadrotor
      Vector3d kp;///< Proportional gains for rpy control
      Vector3d kd;///< Derivative gains for rpy control
      Vector3d a0;///< Offsets in global acceleration due to wind etc
      Vector3d tau0;///< offsets in torque in body frame
      Vector3d acc;///< Internal acceleration computed from last Step function
      //, kp(5,5,5)
      //, kd(7,7,7)
    public:
      QRotorIDModel():System<QRotorIDState,15,4,13>(QRotorIDManifold::Instance())
                      , so3(SO3::Instance())
                      , kt(0.16)
                      , kp(10,10,1.0)
                      , kd(5,5,5)
                      , e3(0,0,1)
                      , a0(0,0,0)
                      , tau0(0,0,0)
                      , g(0,0,-9.81)
                      , acc(0,0,0)
      {
      }
      double Step(QRotorIDState &xb, double t, const QRotorIDState &xa, const Vectorcd &u, double h, const Vectormd *p, Matrixnd *A, Matrixncd *B, Matrixnmd *C)
      {
          Vector3d rpy;
          so3.g2q(rpy,xa.R);
          Vector3d rpy_cmd = xa.u + u.tail<3>()*h;
          for(int i = 0; i < 3; i++)
          {
            rpy_cmd[i] = rpy_cmd[i]>M_PI?(rpy_cmd[i]-2*M_PI):(rpy_cmd[i]<-M_PI)?(rpy_cmd[i]+2*M_PI):rpy_cmd[i];
          }
          Vector3d erpy = rpy - rpy_cmd;
          for(int i = 0; i < 3; i++)
          {
            erpy[i] = erpy[i]>M_PI?(-erpy[i]+2*M_PI):(erpy[i]<-M_PI)?(-erpy[i]-2*M_PI):erpy[i];
          }
          Matrix3d Mrpy;//Convert omega to rpydot
          Mrpy<<1, sin(rpy(0))*tan(rpy(1)), cos(rpy(0))*tan(rpy(1)),
                0, cos(rpy(0)),             -sin(rpy(0)),
                0, sin(rpy(0))*(1/cos(rpy(1))), cos(rpy(0))*(1/cos(rpy(1)));
          Vector3d erpydot = Mrpy*xa.w - u.tail<3>();
          if(p == 0)//No Parameters provided
          {
              acc = (kt*u(0))*(xa.R * e3) + g + xa.R*a0;
              //Translational Part
              xb.v = xa.v + h*(acc);
              //Rotational Part
              Vector3d omegadot = -kp.cwiseProduct(erpy) - kd.cwiseProduct(erpydot) + tau0;
              xb.w = xa.w + omegadot*h;
          }
          else
          {
            //Params p is 13x1 vector with kt, kp, kd, a0, tau0
            acc = ((*p)(0)*u(0))*(xa.R * e3) + g + xa.R*p->segment<3>(7);
            //Translational Part
            xb.v = xa.v + h*(acc);
            //Rotational Part
            Vector3d omegadot = -p->segment<3>(1).cwiseProduct(erpy) - p->segment<3>(4).cwiseProduct(erpydot) + p->tail<3>();
            xb.w = xa.w + omegadot*h;
          }
          xb.p = xa.p + 0.5*(xb.v + xa.v)*h;
          Vector3d w_avg = 0.5*(xb.w + xa.w);
          Matrix3d dR;
          so3.exp(dR,w_avg*h);
          xb.R = xa.R*dR;
          xb.u = rpy_cmd;
          //Specify A,B,C also TODO
      }
  };
}

#endif // QROTORIDMODEL_H
