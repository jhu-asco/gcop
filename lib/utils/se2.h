#ifndef GCOP_SE2_H
#define GCOP_SE2_H

#include "group.h"

namespace gcop {

  using namespace Eigen;

  class SE2 : public Group<3,3> {
  public:

    SE2();

    static SE2& Instance();
  
    void inv(Matrix3d &mi, const Matrix3d &m) const;

    void hat(Matrix3d &vh, const Vector3d &v) const;

    void hatinv(Vector3d &v, const Matrix3d &vh) const;

    void Tg(Matrix3d& m, const Matrix3d &g) const;
    
    void Ad(Matrix3d& m, const Matrix3d &g) const;

    void ad(Matrix3d &m, const Vector3d &v) const;

    void adinv(Vector3d& v, const Matrix3d& m) const;

    void exp(Matrix3d &m, const Vector3d &v) const;

    void log(Vector3d& v, const Matrix3d& m) const;

    void plog(Vector3d& v, const Matrix3d& m) const;

    void cay(Matrix3d& g, const Vector3d &v) const;
    
    void cayinv(Vector3d& v, const Matrix3d& m) const;

    void dcay(Matrix3d& m, const Vector3d& v) const;
    
    void dcayinv(Matrix3d& m, const Vector3d& v) const;

    void q2g(Matrix3d &g, const Vector3d &q) const;

    void g2q(Vector3d &q, const Matrix3d &g) const;

    double tol;               ///< numerical tolerance
    
  };
}

#endif
