#ifndef GCOP_IMUMANIFOLD_H
#define GCOP_IMUMANIFOLD_H

#include "manifold.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 9, 1> Vector9d;
  typedef Matrix<double, 9, 9> Matrix9d;

  struct ImuState {
  ImuState() : 
    R(Matrix3d::Identity()), 
      bg(Vector3d::Zero()), 
      ba(Vector3d::Zero()), 
      P(Matrix9d::Identity()) {
    }
    
    Matrix3d R;   ///< rotation matrix
    Vector3d bg;  ///< gyro bias
    Vector3d ba;  ///< acceleration bias
    
    Matrix9d P;   ///< covariance
  };
  
  //  typedef pair<Matrix3d, Vector6d> ImuState;
  
  class ImuManifold : public Manifold<ImuState, 9> {
    
  public:
    static ImuManifold& Instance();

    void Lift(Vector9d &v,
              const ImuState &xa,
              const ImuState &xb);      

    void Retract(ImuState &xb, 
                 const ImuState &xa,
                 const Vector9d &v);

    void dtau(Matrix9d &M, const Vector9d &v);

    void Adtau(Matrix9d &M, const Vector9d &v);

  private:
    ImuManifold();
  };  
}


#endif
