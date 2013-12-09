#ifndef GCOP_BODY3D_H
#define GCOP_BODY3D_H

#include "system.h"
#include "body3dmanifold.h"
#include "so3.h"
#include <limits>
#include <iostream>
#include <utility>
#include "function.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 12, 12> Matrix12d;



class Geometry
{
public:
  enum {SPHERE, BOX, CYLINDER, MESH} type;
};

class Sphere : public Geometry
{
public:
  Sphere() { this->clear(); };
  double radius;
  void clear()
  {
    radius = 0;
  };
};

class Box : public Geometry
{
public:
  Box() { this->clear(); };
  Vector3d dim;

  void clear()
  {
    dim << 1, 1, 1;
  };
};

class Cylinder : public Geometry
{
public:
  Cylinder() { this->clear(); };
  double length;
  double radius;

  void clear()
  {
    length = 0;
    radius = 0;
  };
};

class Mesh : public Geometry
{
public:
  Mesh() { this->clear(); };
  std::string filename;
  Vector3d scale;

  void clear()
  {
    filename.clear();
    scale << 1,1,1;
  };
  //  bool initXml(TiXmlElement *);
  //  bool fileExists(std::string filename);
};


class Material
{
public:
  Material() { this->clear(); };
  std::string name;
  std::string texture_filename;
  Vector4d color;

  void clear()
  {
    color << .5, .5, .5, .5;
    texture_filename.clear();
    name.clear();
  };
  //  bool initXml(TiXmlElement* config);
};
/*
class Visual
{
public:
  Visual() { this->clear(); };
  Matrix4d origin;
  boost::shared_ptr<Geometry> geometry;

  std::string material_name;
  boost::shared_ptr<Material> material;

  void clear()
  {
    origin.setIdentity();
    material_name.clear();
    material.reset();
    geometry.reset();
    this->group_name.clear();
  };
  //  bool initXml(TiXmlElement* config);
  std::string group_name;
};


class Collision
{
public:
  Collision() { this->clear(); };
  Matrix4d origin;
  //  Pose origin;
  boost::shared_ptr<Geometry> geometry;

  void clear()
  {
    origin.setIdentity();
    geometry.reset();
    this->group_name.clear();
  };
  //  bool initXml(TiXmlElement* config);
  std::string group_name;
};
*/
  
  /**
   * A single rigid body system
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <int c = 6> class Body3d : public System<Body3dState, Matrix<double, c, 1>, 12, c> {
    
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, c, c> Matrixcd;
    typedef Matrix<double, 6, c> Matrix6xcd;
    typedef Matrix<double, 12, c> Matrix12xcd;
    
  public:
    
    Body3d();
    
    Body3d(const Vector3d &ds, double m);
    
    Body3d(const Vector3d &ds, double m, const Vector3d &J);
    
    
    double Step(Body3dState &xb, double t, const Body3dState &xa, 
                const Vectorcd &u, double h,
                Matrix12d *A = 0, Matrix12xcd *B = 0);

    double Fsymp(double t, Body3dState &xb, const Body3dState &xa, 
                 const Vectorcd &u, double h,
                 Matrix12d *A = 0, Matrix12xcd *B = 0);
    
    double Feuler(double t, Body3dState &xb, const Body3dState &xa, 
                  const Vectorcd &u, double h,
                  Matrix12d *A = 0, Matrix12xcd *B = 0);    
    
    void ID(Vector6d &f,
            double t, const Body3dState &xb, const Body3dState &xa,
            const Vectorcd &u, double h);
    
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

    /// TODO: visual element
    //    boost::shared_ptr<Visual> visual;
    
    /// TODO: collision element
    //    boost::shared_ptr<Collision> collision;
    
    Vector3d ds;    ///< body dimensions
    
    double m;      ///< mass
    
    Vector3d J;    ///< 3x1 inertia matrix
    
    Vector6d I;    ///< 6x1 spatial inertia matrix
    
    Vector3d Dw;   ///< linear rotational damping terms
    
    Vector3d Dv;   ///< linear translational damping terms
    
    Vector3d fp;   ///< constant position force in spatial frame (e.g. due to gravity)
    
    Matrix6xcd Bu; ///< control input transformation

    bool symp;     ///< symplectic?

  };  
  
  
  template <int c> 
    Body3d<c>::Body3d() : 
    System<Body3dState, Matrix<double, c, 1>, 12, c>(Body3dManifold::Instance(), Rn<c>::Instance()),
    ds(1, .5, .5), m(1),
    Dw(0,0,0),
    Dv(0,0,0),
    fp(0,0,0),
    Bu(Matrix<double, 6, c>::Identity()),
    symp(true) {
    
    Compute(J, m, ds);
    Compute(I, m, ds);
  }
  
  
  template <int c> 
    Body3d<c>::Body3d(const Vector3d &ds, double m) : 
    System<Body3dState, Matrix<double, c, 1>, 12, c>(Body3dManifold::Instance(), Rn<c>::Instance()),
    ds(ds), m(m),
    Dw(0,0,0),
    Dv(0,0,0),
    fp(0,0,0),
    Bu(Matrix<double, 6, c>::Identity()),
    symp(true) {

    Compute(J, m, ds);
    Compute(I, m, ds);
  }

  template <int c> 
    Body3d<c>::Body3d(const Vector3d &ds, double m, const Vector3d &J) : 
    System<Body3dState, Matrix<double, c, 1>, 12, c>(Body3dManifold::Instance(), Rn<c>::Instance()),
    ds(ds), m(m), J(J),
    Dw(0,0,0),
    Dv(0,0,0),
    fp(0,0,0),
    Bu(Matrix<double, 6, c>::Identity()),
    symp(true) {
    
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

  
  
  static void Gcay(Matrix3d& m, const Vector3d &w, const Vector3d &J, double h, bool plus = true) {
    SO3 &so3 = SO3::Instance();
    Matrix3d wh;
    so3.hat(wh, w);
    
    Matrix3d Jwh;
    Vector3d Jw = J.cwiseProduct(w);
    so3.hat(Jwh, Jw);
    
    int s = (plus ? 1 : -1);
    
    m =  (so3.Id + ((h*h/2)*w)*w.transpose())*J.asDiagonal() + (s*h/2)*(-Jwh + wh*J.asDiagonal()) + (h*h/4*w.dot(Jw))*so3.Id;
    /*
      m << J1*((h^2*v1^2)/4 + 1) + (J1*h^2*v1^2)/2 + (J2*h^2*v2^2)/4 + (J3*h^2*v3^2)/4,  J2*((v1*v2*h^2)/4 + (v3*h)/2) - (J3*h*v3)/2 + (J2*h^2*v1*v2)/4,              (J2*h*v2)/2 - J3*((h*v2)/2 - (h^2*v1*v3)/4) + (J3*h^2*v1*v3)/4]
      [              (J3*h*v3)/2 - J1*((h*v3)/2 - (h^2*v1*v2)/4) + (J1*h^2*v1*v2)/4, J2*((h^2*v2^2)/4 + 1) + (J1*h^2*v1^2)/4 + (J2*h^2*v2^2)/2 + (J3*h^2*v3^2)/4,              J3*((v2*v3*h^2)/4 + (v1*h)/2) - (J1*h*v1)/2 + (J3*h^2*v2*v3)/4]
      [              J1*((v1*v3*h^2)/4 + (v2*h)/2) - (J2*h*v2)/2 + (J1*h^2*v1*v3)/4,              (J1*h*v1)/2 - J2*((h*v1)/2 - (h^2*v2*v3)/4) + (J2*h^2*v2*v3)/4, J3*((h^2*v3^2)/4 + 1) + (J1*h^2*v1^2)/4 + (J2*h^2*v2^2)/4 + (J3*h^2*v3^2)/2]
      
    */
  }

  template <int c>
    double Body3d<c>::Step(Body3dState &xb, double t, const Body3dState &xa, 
                           const Matrix<double, c, 1> &u, double h,
                           Matrix12d *A, Matrix<double, 12, c> *B) {
    if (symp)
      return Fsymp(t, xb, xa, u, h, A, B);
    else
      return Feuler(t, xb, xa, u, h, A, B);      
  }
  
  template <int c>
    double Body3d<c>::Fsymp(double t, Body3dState &xb, const Body3dState &xa, 
                            const Matrix<double, c, 1> &u, double h,
                            Matrix12d *A, Matrix<double, 12, c> *B) {
    // cwiseProduct
    
    SO3 &so3 = SO3::Instance();
    
    const Matrix3d &Ra = xa.first;  
    const Vector3d &pa = xa.second.head<3>();
    const Vector3d &wa = xa.second.segment<3>(3);
    const Vector3d &va = xa.second.tail<3>();
    
    Matrix3d &Rb = xb.first;
    
    Matrix3d Da, Db;
    
    // initialize wb
    Vector3d wb = wa;
    so3.dcayinv(Da, -h*wa);
    so3.dcayinv(Db, h*wb);
    
    Vector3d Jwa = J.cwiseProduct(wa);
    Vector3d Jwb = J.cwiseProduct(wb);
    
    Vector3d fuw = Bu.topRows(3)*u;
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
    
    Vector3d vb = va + h/m*(fp + Ra*fuv);
    
    Vector3d pb = pa + h*vb;
    
    xb.second.head<3>() = pb;
    xb.second.segment<3>(3) = wb;
    xb.second.tail<3>() = vb;

    //    return 1;
    
    const Matrix3d &I3 = Matrix3d::Identity();
    
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
    double Body3d<c>::Feuler(double t, Body3dState &xb, const Body3dState &xa, 
                             const Matrix<double, c, 1> &u, double h,
                             Matrix12d *A, Matrix<double, 12, c> *B) {
    // cwiseProduct
    
    SO3 &so3 = SO3::Instance();
    
    const Matrix3d &Ra = xa.first;  
    const Vector3d &pa = xa.second.head<3>();
    const Vector3d &wa = xa.second.segment<3>(3);
    const Vector3d &va = xa.second.tail<3>();
    
    Matrix3d &Rb = xb.first;

    Vector3d Jwa = J.cwiseProduct(wa);

    Vector3d wb = wa + h*(Jwa.cross(wa) + Bu.topRows(3)*u).cwiseQuotient(J);
    
    Matrix3d chwb;
    so3.cay(chwb, h*wb);
    Rb = Ra*chwb;
    
    Vector3d fuv = Bu.bottomRows(3)*u;
    
    Vector3d vb = va + h/m*(fp + Ra*Bu.bottomRows(3)*u);
    
    xb.second.head<3>() = pa + h*vb;
    xb.second.segment<3>(3) = wb;
    xb.second.tail<3>() = vb;

    Matrix3d Jwah;
    so3.hat(Jwah, Jwa);

    Matrix3d wah;
    so3.hat(wah, wa);
    
    const Matrix3d &I3 = Matrix3d::Identity();
    
    if (A) {
      Matrix6d D1;
      D1.setZero();
      D1.topLeftCorner<3,3>() = -Matrix3d(J.asDiagonal()) - h*Jwah + wah*(J.asDiagonal());
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

}
#endif
