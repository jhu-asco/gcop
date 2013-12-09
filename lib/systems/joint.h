#ifndef GCOP_JOINT_H
#define GCOP_JOINT_H


namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;

  class Joint {
  public:
    //    Joint();

    Vector6d a;    ///< axis
    Matrix4d gp;   ///< fixed xform from parent to joint
    Matrix4d gc;   ///< fixed xfom from child to joint

    Matrix6d Ac;   ///< change of frame Ad(gc)
    Matrix4d gpi;  ///< from parent to joint inverse
    Matrix4d gci;  ///< from child to joint inverse

    Vector6d S;    ///< Jacobian S=Ac*a
  };
}

#endif
