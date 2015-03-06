#ifndef GCOP_KINBODY3D_H
#define GCOP_KINBODY3D_H

#include "system.h"
#include <utility>
#include "kinbody3dmanifold.h"
#include "rn.h"
#include "se3.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
   /**
   * A simple kinematic control system modeled as a rigid body in 3d. It is meant
   * for parameter estimation problems since there is no dynamics. 
   *
   * Author: Matt Sheckells mshecke1(at)jhu.edu
   */
  template<int _nu = 6> class Kinbody3d : public System<Matrix4d, 6, _nu> {
  public:
    typedef Matrix<double, _nu, 1> Vectorud;
    typedef Matrix<double, 6, _nu> Matrix6ud;
    typedef Matrix<double, 6, 6> Matrix6d;

    Kinbody3d();
    
    virtual double Step(Matrix4d &xb, double t, const Matrix4d &xa,
                const Vectorud &u, double h,  const VectorXd *p = 0,
                Matrix6d *A = 0, Matrix6ud *B = 0, Matrix<double, 6, Dynamic> *C = 0);
    
    
    Vector3d d; ///< body dimensions
    Matrix6ud Bu; 
  };

template <int _nu> 
Kinbody3d<_nu>::Kinbody3d() : System<Matrix4d, 6, _nu>(Kinbody3dManifold::Instance()), 
  d(.1, .05, .05),
  Bu(Matrix<double, 6, _nu>::Identity())
{
}

template <int _nu>
double Kinbody3d<_nu>::Step(Matrix4d &xb, double t, const Matrix4d &xa,
                       const Vectorud &u, double h,  const VectorXd *p,
                       Matrix6d *A, Matrix6ud *B, Matrix<double, 6, Dynamic> *C )
{
  SE3 &se3 = SE3::Instance();
  Matrix4d m;
  Vector6d u_kin = Bu*u;  


  se3.cay(m, h*u_kin);
  xb = xa*m;
  
  if (A) {
    se3.cay(m, -h*u_kin);
    se3.Ad(*A, m);
  }
  
  if (B) {
    Matrix6d mb;
    se3.dcay(mb, -h*u_kin);
    *B = h*mb*Bu;
  }

  if (C)
    C->setZero();
}

} // gcop
#endif
