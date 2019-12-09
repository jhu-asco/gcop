#ifndef GCOP_YAWVELFOLLOW_H
#define GCOP_YAWVELFOLLOW_H

#include "constraint.h"
#include <math.h>
#include "so3.h"

namespace gcop {
  
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic>
    class YawVelFollow : public Constraint<T, _nx, _nu, _np> {
  public:

  typedef Matrix<double, Dynamic, 1> Vectorgd;
  typedef Matrix<double, Dynamic, _nx> Matrixgxd;
  typedef Matrix<double, Dynamic, _nu> Matrixgud;
  typedef Matrix<double, Dynamic, _np> Matrixgpd;

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;


  YawVelFollow();
  
  bool operator()(Vectorgd &g,
                  double t, const T &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrixgxd *dgdx = 0, Matrixgud *dgdu = 0,
                  Matrixgpd *dgdp = 0);    
  
  };
  
  
  template <typename T, int _nx, int _nu, int _np> 
    YawVelFollow<T, _nx, _nu, _np>::YawVelFollow() :
    Constraint<T, _nx, _nu, _np>(1)
  {
    this->ng = 1;
  }
  
  template <typename T, int _nx, int _nu, int _np> 
    bool YawVelFollow<T, _nx, _nu, _np>::operator()(Vectorgd &g,
                                             double t, const T &x, const Vectorcd &u,
                                             const Vectormd *rho, 
                                             Matrixgxd *dgdx, Matrixgud *dgdu,
                                             Matrixgpd *dgdp)
    {
      g.resize(1);
      const Vector3d &vel = x.v; // velocity
      const Matrix3d &Rot = x.R; // Rotation Matrix
     
      double yaw = atan2(Rot(1,0),Rot(0,0));//Yaw Based on SO3 Code
      //NOTE: If both arguments are 0, yaw is undefined, but atan2 returns 0
      double yaw_d = atan2(vel[1],vel[0]);//Velocity following yaw
      
      double vsq = vel[1]*vel[1] + vel[0]*vel[0]; //Velocity magnitude
      double d;
      if (vsq < 1e-10) { //At small velocity, constraint is ill-defined
        d = 0;
      } else {
        double d = yaw - yaw_d;  // distance between angles
        if (d > M_PI){
          d -= 2*M_PI;
        }
        if (d < -M_PI){
          d += 2*M_PI;
        }
      }

      g[0] = d*d; // must be 0 for non-collision
      if (dgdx) {
        dgdx->resize(1, _nx);
        dgdx->setZero();
        //THIS IS A BUG.  For states that are not Body3dState, it won't necessarily work.
        //It assumes that the position indices of the x vector are 3,4,5, but other states use 0,1,2
        //DO NOT TRUST 
        if ((abs(d) < 1e-10)||(vsq < 1e-10)) {//Derivative is 0
          return true;
        }
        double R1 = Rot(1,0);
        double R0 = Rot(0,0);
        double R_sq = R0*R0 + R1*R1;
        if (R_sq < 1e-10) {//Derivative is 0
          return true;
        }
        Vector3d eR;
        SO3::Instance().log(eR,Rot);
        double theta = eR.norm();
        if (theta < 1e-10) {//Derivative is 0
          return true;
        }
        double partial_dx = 2*d;
        Vector3d dThetadeR = eR/theta;
        Vector3d temp(0,2*eR[1],2*eR[2]);
        Vector3d dR0deR = (eR[1]*eR[1] + eR[2]*eR[2])*(-sin(theta))*dThetadeR + (cos(theta)-1)*temp;
        Vector3d temp2(-eR[1]*(cos(theta)-1),-eR[0]*(cos(theta)-1),sin(theta));
        Vector3d dR1deR = eR[2]*cos(theta)*dThetadeR + temp2 + eR[0]*eR[1]*sin(theta)*dThetadeR;
        //0,1,2 are log(R), 9,10,11 are velocity
        dgdx->row(0).head(3) = partial_dx*((R0/R_sq)*dR1deR - (R1/R_sq)*dR0deR);
        dgdx->row(0)[9] = partial_dx*(vel[1]/vsq);
        dgdx->row(0)[10] = partial_dx*(-vel[0]/vsq);
      }
    }
};



#endif
