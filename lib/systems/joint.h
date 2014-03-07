#ifndef GCOP_JOINT_H
#define GCOP_JOINT_H

#include <Eigen/Dense>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
  
  class Joint {
  public:
    Joint();
    
    Joint(const Vector6d &a,
          const Matrix4d &gp,
          const Matrix4d &gc,
          double damping = 0,
          double friction = 0);
    
    virtual ~Joint();
    
    void Init();

    Vector6d a;    ///< axis
    Matrix4d gp;   ///< fixed xform from parent to joint
    Matrix4d gc;   ///< fixed xfom from child to joint

    double damping;   ///< damping terms: Ns/m for trans or Nm s/rad for rev
    double friction;  ///< static friction: N for prismatic, Nm for rev

    Matrix6d Ac;   ///< change of frame Ad(gc)
    Matrix4d gpi;  ///< from parent to joint inverse
    Matrix4d gci;  ///< from child to joint inverse

    Vector6d S;    ///< Jacobian S=Ac*a
		string name;   ///< Unique name of the joint
  };
}

#endif
