#ifndef GCOP_INSMANIFOLD_H
#define GCOP_INSMANIFOLD_H

#include "manifold.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 15, 1> Vector15d;
  typedef Matrix<double, 15, 15> Matrix15d;

  struct InsState {
  InsState() : 
    R(Matrix3d::Identity()), 
      bg(Vector3d::Zero()), 
      ba(Vector3d::Zero()), 
      p(Vector3d::Zero()),
      v(Vector3d::Zero()),
      P(Matrix15d::Identity()) {
    }
    
    Matrix3d R;   ///< rotation matrix
    Vector3d bg;  ///< gyro bias
    Vector3d ba;  ///< acceleration bias    

    Vector3d p;   ///< position
    Vector3d v;   ///< velocity
    
    Matrix15d P;   ///< covariance
  };
  
  //  typedef pair<Matrix3d, Vector6d> InsState;
  
  class InsManifold : public Manifold<InsState, 15> {
    
  public:
    static InsManifold& Instance();

    void Lift(Vector15d &v,
              const InsState &xa,
              const InsState &xb);      

    void Retract(InsState &xb, 
                 const InsState &xa,
                 const Vector15d &v);

    void dtau(Matrix15d &M, const Vector15d &v);

    void Adtau(Matrix15d &M, const Vector15d &v);

  private:
    InsManifold();
  };  
}


#endif
