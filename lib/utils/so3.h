#ifndef GCOP_SO3_H
#define GCOP_SO3_H

#include "group.h"

namespace gcop {

  using namespace Eigen;

  class SO3 : public Group<3,3> {
  public:

    SO3();

    static SO3& Instance();
  
    void inv(Matrix3d &mi, const Matrix3d &m) const;

    void hat(Matrix3d &vh, const Vector3d &v) const;

    void hatinv(Vector3d &v, const Matrix3d &vh) const;

    void Tg(Matrix3d& m, const Matrix3d &g) const;
    
    void Ad(Matrix3d& m, const Matrix3d &g) const;

    void ad(Matrix3d &m, const Vector3d &v) const;

    void adinv(Vector3d& v, const Matrix3d& m) const;

    void exp(Matrix3d &m, const Vector3d &v) const;

    void log(Vector3d& v, const Matrix3d& m) const;

    void cay(Matrix3d& g, const Vector3d &v) const;
    
    void cayinv(Vector3d& v, const Matrix3d& m) const;

    void dcay(Matrix3d& m, const Vector3d& v) const;
    
    void dcayinv(Matrix3d& m, const Vector3d& v) const;

    void dexp(Matrix3d& m, const Vector3d& v) const;
    
    void dexpinv(Matrix3d& m, const Vector3d& v) const;

    void skew(Vector3d& v, const Matrix3d& m) const;  

    double roll(const Matrix3d& m) const;

    double pitch(const Matrix3d& m) const;

    double yaw(const Matrix3d& m) const;

    void q2g(Matrix3d &m, const Vector3d &rpy);

    void quat2g(Matrix3d &m, const Vector4d &wxyz);

    void g2quat(Vector4d &wxyz, const Matrix3d &m);

    void g2q(Vector3d &rpy, const Matrix3d &m);

    double tol;               ///< numerical tolerance
    
  };
}

#endif
