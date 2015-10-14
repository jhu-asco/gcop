#ifndef GCOP_POSEMANIFOLD_H
#define GCOP_POSEMANIFOLD_H

#include "manifold.h"
#include <iostream>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;

  struct PoseState {
  PoseState() : 
    R(Matrix3d::Identity()), 
      p(Vector3d::Zero()),
      P(Matrix6d::Identity()) {
    }
    
    Matrix3d R;   ///< rotation matrix
    Vector3d p;   ///< position
    
    Matrix6d P;   ///< covariance
  void resize(int n){
      std::cout << "Warning: type is PoseState, resize is NOT implemented" << std::endl;
    }
  };
  
  //  typedef pair<Matrix3d, Vector6d> InsState;
  
  class PoseManifold : public Manifold<PoseState, 6> {
    
  public:
    static PoseManifold& Instance();

    void Lift(Vector6d &v,
              const PoseState &xa,
              const PoseState &xb);      

    void Retract(PoseState &xb, 
                 const PoseState &xa,
                 const Vector6d &v);

    void dtau(Matrix6d &M, const Vector6d &v);

    void Adtau(Matrix6d &M, const Vector6d &v);

  private:
    PoseManifold();
  };  
}


#endif
