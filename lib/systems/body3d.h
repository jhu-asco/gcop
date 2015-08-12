#ifndef GCOP_BODY3D_H
#define GCOP_BODY3D_H

#include "system.h"
#include "body3dmanifold.h"
#include "so3.h"
#include "se3.h"
#include <limits>
#include <iostream>
#include <utility>
#include "function.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 12, 1> Vector12d;
  typedef Matrix<double, 12, 12> Matrix12d;
  typedef Matrix<double, 12, 6> Matrix12x6d;
  /**
   * A single rigid body system
   *
   * Note: damping is currently supported only by Euler's method
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <int c = 6> class Body3d : public System<Body3dState, 12, c> {
  protected:
  
  typedef Matrix<double, c, 1> Vectorcd;
  typedef Matrix<double, c, c> Matrixcd;
  typedef Matrix<double, 6, c> Matrix6xcd;
  typedef Matrix<double, 12, c> Matrix12xcd;
  typedef Matrix<double, 12, Dynamic> Matrixnmd;

  public
  
  Body3d(int np =0);
  
  Body3d(const Vector3d &ds, double m, int np = 0);
  
  Body3d(const Vector3d &ds, double m, const Vector3d &J, int np = 0);
  
  virtual ~Body3d();
  
  virtual double Step(Body3dState &xb, double t, const Body3dState &xa, 
                      const Vectorcd &u, double h, const VectorXd *p = 0,
                      Matrix12d *A = 0, Matrix12xcd *B = 0, Matrixnmd *C = 0);

  virtual bool NoiseMatrix(Matrix12d &L, double t, const Body3dState &x, 
                           const Vectorcd &u, double h, const VectorXd *p = 0);

  virtual bool Noise(Matrix12d &Q, double t, const Body3dState &x, 
                     const Vectorcd &u, double h, const VectorXd *p = 0);
  
  
  double SympStep(double t, Body3dState &xb, const Body3dState &xa, 
                  const Vectorcd &u, double h, const VectorXd *p = 0,
                  Matrix12d *A = 0, Matrix12xcd *B = 0, Matrixnmd *C = 0);
  
  double EulerStep(double t, Body3dState &xb, const Body3dState &xa, 
                   const Vectorcd &u, double h, const VectorXd *p = 0,
                   Matrix12d *A = 0, Matrix12xcd *B = 0, Matrixnmd *C = 0);    
  
  

  /**
   * Inverse dynamics
   * @param f inverse dynamics vector
   * @param t time
   * @param xb end state
   * @param xa start state
   * @param u control
   * @param h time-step
   */
  void ID(Vector6d &f,
          double t, const Body3dState &xb, const Body3dState &xa,
          const Vectorcd &u, double h);
  
  void StateAndControlsToFlat(VectorXd &y, const Body3dState &x, const Vectorcd &u);
  
  void FlatToStateAndControls(Body3dState &x, Vectorcd &u, const std::vector<VectorXd> &y);
  
  /**
   * Compute moments of inertia J given mass m and dimensions ds
   * @param J moments of inertia (3x1 vector)
   * @param m mass
   * @param ds dimensions (assume a rectangular body)
   */
  static void Compute(Vector3d &J, double m, const Vector3d &ds);
  
  /**
   * Compute moments of inertia tensor I given mass m and dimensions ds
   * @param J moments of inertia (6x1 vector)
   * @param m mass
   * @param ds dimensions (assume a rectangular body)
   */
  static void Compute(Vector6d &I, double m, const Vector3d &ds);

  Vector3d ds;    ///< body dimensions
  
  double m;      ///< mass
  
  Vector3d J;    ///< 3x1 principal moments of inertia
  
  Vector6d I;    ///< 6x1 spatial inertia matrix

  Vector6d Ii;   ///< 6x1 spatial inertia matrix inverse
  
  Vector3d Dw;   ///< linear rotational damping terms
  
  Vector3d Dv;   ///< linear translational damping terms
    
  Vector3d fp;   ///< constant position force in spatial frame (e.g. due to gravity)
  
  Matrix6xcd Bu; ///< control input transformation
  
  bool symp;     ///< symplectic?
  
  string name;   ///< Unique name of the body

  double sw;     ///< process noise: standard deviation in torque (zero by default)
  double sv;     ///< process noise: standard deviation in force (zero by default)


  bool constantVelocity;  ///< whether to ignore the rigid body dynamics, and use a simple constant velocity model, in which the inputs are ignored, and it is assumed that the body moves with constant velocity subject to a white noise acceleration (false by default)
  };  
  


  template <int c> 
    Body3d<c>::Body3d(int np) : 
    System<Body3dState, 12, c>(Body3dManifold::Instance(), c, np),
    ds(.6, .6, .2), m(1),
    Dw(0,0,0),
    Dv(0,0,0),
    fp(0,0,0),
    Bu(Matrix<double, 6, c>::Identity()),
    symp(true),
    sw(0), sv(0), constantVelocity(false) {
    
    Compute(J, m, ds);
    Compute(I, m, ds);
  }
  
  
  template <int c> 
    Body3d<c>::Body3d(const Vector3d &ds, double m, int np) : 
    System<Body3dState, 12, c>(Body3dManifold::Instance(),c, np),
    ds(ds), m(m),
    Dw(0,0,0),
    Dv(0,0,0),
    fp(0,0,0),
    Bu(Matrix<double, 6, c>::Identity()),
    symp(true),
    sw(0), sv(0), constantVelocity(false) {

    Compute(J, m, ds);
    Compute(I, m, ds);
  }

  template <int c> 
    Body3d<c>::Body3d(const Vector3d &ds, double m, const Vector3d &J, int np) : 
    System<Body3dState, 12, c>(Body3dManifold::Instance(), c, np),
    ds(ds), m(m), J(J),
    Dw(0,0,0),
    Dv(0,0,0),
    fp(0,0,0),
    Bu(Matrix<double, 6, c>::Identity()),
    symp(true),
    sw(0), sv(0), constantVelocity(false) {
    
  }

template<int c>  Body3d<c>::~Body3d()
{
}

  template <int c> 
    void Body3d<c>::Compute(Vector3d &J, double m, const Vector3d &ds) {
    J[0] = m*(ds[1]*ds[1] + ds[2]*ds[2])/3;
    J[1] = m*(ds[0]*ds[0] + ds[2]*ds[2])/3;
    J[2] = m*(ds[0]*ds[0] + ds[1]*ds[1])/3;
  }

  template <int c> 
    void Body3d<c>::Compute(Vector6d &I, double m, const Vector3d &ds) {
    I[0] = m*(ds[1]*ds[1] + ds[2]*ds[2])/3;
    I[1] = m*(ds[0]*ds[0] + ds[2]*ds[2])/3;
    I[2] = m*(ds[0]*ds[0] + ds[1]*ds[1])/3;
    I[3] = m; 
    I[4] = m;
    I[5] = m;
  }


  template <int c> 
    bool Body3d<c>::NoiseMatrix(Matrix12d &L, 
                                double t, const Body3dState &x, 
                                const Vectorcd &u, double dt, 
                                const VectorXd *p) {    
    // @mk: you could use a slighly more accurate noise
    L.setZero();
    double dt2 = dt*dt/2;
    L(0,6) = sw*dt2/J(0);
    L(1,7) = sw*dt2/J(1);
    L(2,8) = sw*dt2/J(2);
    L(3,9) = sv*dt2/m;
    L(4,10) = sv*dt2/m;
    L(5,11) = sv*dt2/m;
    
    L(6,6) = sw*dt/J(0);
    L(7,7) = sw*dt/J(1);
    L(8,8) = sw*dt/J(2);
    L(9,9) = sv*dt/m;
    L(10,10) = sv*dt/m;
    L(11,11) = sv*dt/m;

    return true;
  }


  template <int c> 
    bool Body3d<c>::Noise(Matrix12d &Q, double t, const Body3dState &x, 
                          const Vectorcd &u, 
                          double dt, const VectorXd *p) {
    
    if (constantVelocity) {
      Q.setZero();
      double dt2 = dt*dt;
      double dt3 = dt2*dt;
      
      // a coarse Euler-type approximation to rotational error propagation
      double sw2 = sw*sw;
      Q.topLeftCorner<3,3>().diagonal().setConstant(sw2*dt3/3);
      Q.block<3,3>(0,3).diagonal().setConstant(sw2*dt2/2);
      Q.block<3,3>(3,0).diagonal().setConstant(sw2*dt2/2);
      Q.block<3,3>(3,3).diagonal().setConstant(sw2*dt);    
      
      // exact error propagation for position
      double sv2 = sv*sv;
      Q.block<3,3>(6,6).diagonal().setConstant(sv2*dt3/3);
      Q.block<3,3>(6,9).diagonal().setConstant(sv2*dt2/2);
      Q.block<3,3>(9,6).diagonal().setConstant(sv2*dt2/2);
      Q.block<3,3>(9,9).diagonal().setConstant(sv2*dt);
    } else {
      
    }

    return true;
  }
  
  
  static void Gcay(Matrix3d& m1, const Vector3d &w, const Vector3d &J, double h, bool plus = true) {
    SO3 &so3 = SO3::Instance();
    Matrix3d wh;
    so3.hat(wh, w);
    
    Matrix3d Jwh;
    Vector3d Jw = J.cwiseProduct(w);
    so3.hat(Jwh, Jw);
    
    int s = (plus ? 1 : -1);
    
    m1 =  (so3.Id + ((h*h/2)*w)*w.transpose())*J.asDiagonal() + (s*h/2)*(-Jwh + wh*J.asDiagonal()) + (h*h/4*w.dot(Jw))*so3.Id;
    /*
      m << J1*((h^2*v1^2)/4 + 1) + (J1*h^2*v1^2)/2 + (J2*h^2*v2^2)/4 + (J3*h^2*v3^2)/4,  J2*((v1*v2*h^2)/4 + (v3*h)/2) - (J3*h*v3)/2 + (J2*h^2*v1*v2)/4,              (J2*h*v2)/2 - J3*((h*v2)/2 - (h^2*v1*v3)/4) + (J3*h^2*v1*v3)/4]
      [              (J3*h*v3)/2 - J1*((h*v3)/2 - (h^2*v1*v2)/4) + (J1*h^2*v1*v2)/4, J2*((h^2*v2^2)/4 + 1) + (J1*h^2*v1^2)/4 + (J2*h^2*v2^2)/2 + (J3*h^2*v3^2)/4,              J3*((v2*v3*h^2)/4 + (v1*h)/2) - (J1*h*v1)/2 + (J3*h^2*v2*v3)/4]
      [              J1*((v1*v3*h^2)/4 + (v2*h)/2) - (J2*h*v2)/2 + (J1*h^2*v1*v3)/4,              (J1*h*v1)/2 - J2*((h*v1)/2 - (h^2*v2*v3)/4) + (J2*h^2*v2*v3)/4, J3*((h^2*v3^2)/4 + 1) + (J1*h^2*v1^2)/4 + (J2*h^2*v2^2)/4 + (J3*h^2*v3^2)/2]
      
    */
  }

  template <int c>
    double Body3d<c>::Step(Body3dState &xb, double t, const Body3dState &xa, 
                           const Matrix<double, c, 1> &u, double dt, const VectorXd *p,
                           Matrix12d *A, Matrix<double, 12, c> *B, Matrixnmd *C) {

    // in a constant velocity model the dynamics is ignored:
    if (constantVelocity) {

      Matrix3d dR;
      SO3::Instance().exp(dR, dt*xa.w);
      xb.R = xa.R*dR;
      xb.p = xa.p + dt*xa.v;
      xb.w = xa.w;
      xb.v = xa.v;
      

      if (A) {
        Matrix3d D;
        SO3::Instance().dexp(D, -dt*xa.w);
        A->setIdentity();
        A->topLeftCorner<3,3>() = dR.transpose();
        A->block<3,3>(0,6) = dt*D;           // dR wrt w (trivialized)        
        A->block<3,3>(3,9).diagonal().setConstant(dt);  // dp drt v        
      }
      if (B) {
        B->setZero();
      }      
      
      return 0;
    }
    
    // for regular rigid body dynamics use either a symplectic or a semi-implicit Euler method
    if (symp)
      return SympStep(t, xb, xa, u, dt, p, A, B, C);
    else
      return EulerStep(t, xb, xa, u, dt, p, A, B, C);      
  }
  
  template <int c>
    double Body3d<c>::SympStep(double t, Body3dState &xb, const Body3dState &xa, 
                               const Matrix<double, c, 1> &u, double h, const VectorXd *p,
                               Matrix12d *A, Matrix<double, 12, c> *B, Matrixnmd *C) {
    // cwiseProduct
    
    SO3 &so3 = SO3::Instance();
    
    const Matrix3d &Ra = xa.R;  
    const Vector3d &pa = xa.p;
    const Vector3d &wa = xa.w;
    const Vector3d &va = xa.v;
    
    Matrix3d &Rb = xb.R;
    
    Matrix3d Da, Db;
    
    // initialize wb
    Vector3d wb = wa;
    so3.dcayinv(Da, -h*wa);
    so3.dcayinv(Db, h*wb);
    
    Vector3d Jwa = J.cwiseProduct(wa);
    Vector3d Jwb = J.cwiseProduct(wb);
    
    Vector3d fuw = Bu.topRows(3)*u;
    //Adding external force through parameters:
     if(p != 0)
    {
      assert(p->size() == 6);
      fuw += (*p).head<3>();
    }
    Vector3d fwd = Db.transpose()*Jwb - Da.transpose()*Jwa - h*fuw;
    
    Matrix3d Gp;
    Gcay(Gp, wb, J, h);
    
    Matrix3d Gm;
    Gcay(Gm, wa, J, h, false);
    
    // Newton step
    Matrix3d Gi = Gp.inverse();
    wb = wb - Gi*fwd;
    //    wb = wa + Jwa.cross(wa) + h*fuw;

    // update Gi
    Gcay(Gp, wb, J, h);
    Gi = Gp.inverse();
    
    Matrix3d chwb;
    so3.cay(chwb, h*wb);
    Rb = Ra*chwb;
    
    Vector3d fuv = Bu.bottomRows(3)*u;
    //Adding external force through parameters:
    /*if(p != 0)
    {
      assert(p->size() == 6);
      fuv += Bu.bottomRows(3)*(*p);
    }
    */

    const Matrix3d &I3 = Matrix3d::Identity();

    Vector3d vb;
    Vector3d pb;

    bool vdamp = (Dv.norm() > 1e-16); // is there linear velocity damping
    // modified mass matrix to include damping
    Matrix3d Mp;
    Matrix3d Mm;    

    if (vdamp) {
      Mp = I3*m + h/2*(Ra*Dv.asDiagonal());
      Mm = I3*m - h/2*(Ra*Dv.asDiagonal());

      if(p != 0)
      {
        vb = Mp.inverse()*(Mm*va + h*(fp + Ra*fuv + (*p).tail<3>()));    
      }
      else
      {
        vb = Mp.inverse()*(Mm*va + h*(fp + Ra*fuv));    
      }
      pb = pa + h*vb;
    } else {
      vb = va + h/m*(fp + Ra*fuv);    
      if(p != 0)
      {
        vb = va + h/m*(fp + Ra*fuv + (*p).tail<3>());    
      }
      else
      {
        vb = va + h/m*(fp + Ra*fuv);    
      }
      pb = pa + h*vb;
    }
    
    xb.p = pb;
    xb.w = wb;
    xb.v = vb;

    //    return 1;
        
    if (A) {
      Matrix6d D1;
      D1.setZero();
      D1.topLeftCorner<3,3>() = -Gm;
      D1.bottomRightCorner<3,3>() = -m*I3;
      
      Matrix6d D2;
      D2.setZero();
      Matrix3d fuvh; 
      so3.hat(fuvh, fuv);
      D2.bottomLeftCorner<3,3>() = h*Ra*fuvh;
      
      Matrix6xcd D3; 
      D3.topRows(3) = -h*Bu.topRows(3);
      D3.bottomRows(3) = -h*Ra*Bu.bottomRows(3);

      Matrix3d cb;
      so3.cay(cb, -h*wb);
      
      Matrix3d Adcb;
      so3.Ad(Adcb, cb);
      
      Matrix3d dcb;
      so3.dcay(dcb, -h*wb);
      
      Matrix12d A1;
      A1.setZero();
      A1.topLeftCorner<3,3>() = Adcb;
      A1.block<3,3>(3,3) = I3;
      A1.block<3,3>(0,6) = h*dcb*Gi;
      A1.block<3,3>(3,9) = (h/m)*I3;
      
      A1.block<3,3>(6,6) = Gi;
      A1.block<3,3>(9,9) = I3/m;
      
      Matrix12d A2;
      A2.setIdentity();
      A2.block<6,6>(6,0) = -D2;
      A2.block<6,6>(6,6) = -D1;
      
      (*A) = A1*A2;
      
      if (B) {
        (*B) = -A1.rightCols<6>()*D3;
      }
    }
    return 1;
  }


  template <int c>
    double Body3d<c>::EulerStep(double t, Body3dState &xb, const Body3dState &xa, 
                                const Matrix<double, c, 1> &u, double h, const VectorXd *p,
                                Matrix12d *A, Matrix<double, 12, c> *B, Matrixnmd *C) {
    // cwiseProduct
    
    SO3 &so3 = SO3::Instance();
    
    const Matrix3d &Ra = xa.R;
    const Vector3d &pa = xa.p;
    const Vector3d &wa = xa.w;
    const Vector3d &va = xa.v;
    
    Matrix3d &Rb = xb.R;

    Vector3d Jwa = J.cwiseProduct(wa);

    Vector3d wb; 
    if(p!= 0)
    {
      wb = wa + h*(Jwa.cross(wa) + Dw.cwiseProduct(wa) + Bu.topRows(3)*u + (*p).tail<3>()).cwiseQuotient(J);//Inside brackets + force
    }
    else
    {
      wb = wa + h*(Jwa.cross(wa) + Dw.cwiseProduct(wa) + Bu.topRows(3)*u).cwiseQuotient(J);//Inside brackets + force
    }
    
    Matrix3d chwb;
    so3.cay(chwb, h*wb);
    Rb = Ra*chwb;
    
    Vector3d fuv = Dv.cwiseProduct(Ra.transpose()*va) + Bu.bottomRows(3)*u;
    /*//Adding external force through parameters:
    if(p != 0)
    {
      assert(p->size() == 6);
      fuv += Bu.bottomRows(3)*(*p);
    }
    */

    Vector3d vb = va + h/m*(fp + Ra*fuv);
    if(p!= 0)
    {
      vb += h/m*((*p).tail<3>());
    }
    
    xb.p = pa + h*vb;
    xb.w = wb;
    xb.v = vb;

    Matrix3d Jwah;
    so3.hat(Jwah, Jwa);

    Matrix3d wah;
    so3.hat(wah, wa);
    
    const Matrix3d &I3 = Matrix3d::Identity();
    
    if (A) {
      Matrix6d D1;
      D1.setZero();
      D1.topLeftCorner<3,3>() = -Matrix3d(J.asDiagonal()) - h*(Jwah + Matrix3d(Dw.asDiagonal())) + wah*(J.asDiagonal());
      D1.bottomRightCorner<3,3>() = -m*I3;
      
      Matrix6d D2;
      D2.setZero();
      Matrix3d fuvh;
      so3.hat(fuvh, fuv);
      D2.bottomLeftCorner<3,3>() = h*Ra*fuvh;
      
      Matrix6xcd D3;
      D3.topRows(3) = -h*Bu.topRows(3);
      D3.bottomRows(3) = -h*Ra*Bu.bottomRows(3);

      Matrix3d cb;
      so3.cay(cb, -h*wb);
      
      Matrix3d Adcb;
      so3.Ad(Adcb, cb);
      
      Matrix3d dcb;
      so3.dcay(dcb, -h*wb);
      
      Matrix12d A1;
      A1.setZero();
      A1.topLeftCorner<3,3>() = Adcb;
      A1.block<3,3>(3,3) = I3;
      A1.block<3,3>(0,6) = h*dcb*J.asDiagonal().inverse();
      A1.block<3,3>(3,9) = (h/m)*I3;
      
      A1.block<3,3>(6,6) = J.asDiagonal().inverse();
      A1.block<3,3>(9,9) = I3/m;
      
      Matrix12d A2;
      A2.setIdentity();
      A2.block<3,3>(0,0).diagonal() -= h*Dw.cwiseQuotient(J);
      A2.block<3,3>(3,3) -= (h/m)*Ra*Dv.asDiagonal()*Ra;
      A2.block<6,6>(6,0) = -D2;
      A2.block<6,6>(6,6) = -D1;
      
      (*A) = A1*A2;
      
      if (B) {
        (*B) = -A1.rightCols<6>()*D3;
      }
    }
    return 1;
  }




  template <int c>   
    void Body3d<c>::ID(Vector6d &f,
                       double t, const Body3dState &xb, const Body3dState &xa,
                       const Vectorcd &u, double h) {
  }

  template <int c>
    void Body3d<c>::StateAndControlsToFlat(VectorXd &y, const Body3dState &x,
               const Vectorcd &u) {
    SO3& so3 = SO3::Instance();

    // Flat outputs are x,y,z, and rpy
    y.resize(6);
    y.head<3>() = x.p;
    y(3) = so3.roll(x.R);
    y(4) = so3.pitch(x.R);
    y(5) = so3.yaw(x.R);
  }

  template <int c>
    void Body3d<c>::FlatToStateAndControls(Body3dState &x, Vectorcd &u,
               const std::vector<VectorXd> &y) {
    assert(y.size() >= 2);

    SO3& so3 = SO3::Instance();

    VectorXd y0 = y[0];
    VectorXd y1 = y[1];

    //    x.second.setZero();
    u.setZero();

    so3.q2g(x.R, y0.tail<3>());

    x.p = y0.head<3>();
    x.v = y1.head<3>();
    // TODO: fill in angular velocity and controls
    x.w[0] = 0; x.w[1] = 0;
    x.w[2] = y1(5);
  }
}
#endif
