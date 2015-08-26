#ifndef GCOP_BODY2DFORCE_H
#define GCOP_BODY2DFORCE_H

#include "force.h"
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef pair<Matrix3d, Vector3d> Body2dState;
  typedef Matrix<double, 3, 6> Matrix36d;

  template<int c = 3>
    class Body2dForce : public Force<Body2dState, 3, 6, c> {
  public:

  typedef Matrix<double, c, 1> Vectorcd;
  typedef Matrix<double, 6, c> Matrix6cd;
  typedef Matrix<double, 3, c> Matrix3cd;

  Body2dForce(bool fextParam = false, bool dampParam = false);
  
  virtual void Set(Vector3d &f, const Body2dState &x, double t, 
                   const Vectorcd &u, double h, const VectorXd *p = 0,
                   Matrix36d *A = 0, Matrix3cd *B = 0, Matrix<double, 3, Dynamic> *C = 0);
  
  Vector3d D;      ///< damping terms
  Matrix3cd B;     ///< constant control input transformation
  Vector3d fext;   ///< constant external force

  bool fextParam;  ///< treat the (x,y) forces in fext as a parameter
  bool dampParam;  
  };

  template<int c>
    Body2dForce<c>::Body2dForce(bool fextParam, bool dampParam) : 
    D(0, 0, 0), fext(0, 0, 0), fextParam(fextParam), dampParam(dampParam)
    {
      B.setIdentity();
    }
  
  
  template<int c>
    void Body2dForce<c>::Set(Vector3d &f, const Body2dState &x, double t, 
                             const Vectorcd &u, double h,  const VectorXd *p,
                             Matrix36d *A, Matrix3cd *B, Matrix<double, 3, Dynamic> *C)
    {
      if (fextParam && p) {
        fext.tail<2>() = *p;
        assert(p->size() == 2);
      }

      if (dampParam && p) {
        D = *p;
      }
      
      Vector3d fb; // body-fixed external force
      fb[0] = fext[0];
      fb.tail<2>() = x.first.topLeftCorner<2,2>().transpose()*fext.tail<2>();    

      const Vector3d& v = x.second;
      f = this->B*u - D.cwiseProduct(v.cwiseAbs().cwiseProduct(v)) + fb;
      
      if (A) {
        /*
        A->leftCols<3>().setZero();    
        (*A)(1,0) = -fb[2]; (*A)(2,0)= fb[1]; //SE2::r2hat(fb.tail<2>()).tanspose();    
        A->rightCols<3>() = Matrix3d((-D).asDiagonal());
        */
      }
      
      if (B) {
        //        *B = this->B;
      }
      
      if (C) {
        /*
        if (fextParam && p) {
          C->topRows<1>().setZero();
          C->bottomRows<2>() = x.first.topLeftCorner<2, 2>().transpose();
        } else {
          // C->setZero();
        }
        */
      }
    }

}


#endif
