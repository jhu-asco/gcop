#include "body2dforce.h"
#include "se2.h"
#include <iostream>
#include <assert.h>
#include <limits>

using namespace gcop;
using namespace Eigen;

Body2dForce::Body2dForce(bool fextParam) : D(0, 0, 0), fext(0, 0, 0), fextParam(fextParam)
{
}


void Body2dForce::Set(Vector3d &f, const M3V3d &x, double t, 
                      const Vector3d &u, double h,  const VectorXd *p,
                      Matrix36d *A, Matrix3d *B, Matrix<double, 3, Dynamic> *C )
{
  if (fextParam && p) {
    fext.tail<2>() = *p;
    assert(p->size() == 2);
  }

  Vector3d fb; // body-fixed external force
  fb[0] = fext[0];
  fb.tail<2>() = x.first.topLeftCorner<2,2>().transpose()*fext.tail<2>();    
  f = u - D.cwiseProduct(x.second) + fb;
  
  if (A) {
    A->leftCols<3>().setZero();    
    (*A)(1,0) = -fb[2]; (*A)(2,0)= fb[1]; //SE2::r2hat(fb.tail<2>()).tanspose();    
    A->rightCols<3>() = Matrix3d((-D).asDiagonal());
  }
  
  if (B) {
    B->setIdentity();
  }

  if (C) {
    if (fextParam && p) {
      C->topRows<1>().setZero();
      C->bottomRows<2>() = x.first.topLeftCorner<2, 2>().transpose();
    } else {
      C->setZero();
    }
  }
}
