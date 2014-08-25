#include <limits>
#include "body2d.h"
#include "body2dmanifold.h"
#include "se2.h"
#include "rn.h"
#include <iostream>
#include <assert.h>
#include "body2dforce.h"


using namespace gcop;
using namespace Eigen;

Body2d::Body2d(Body2dForce *force) : 
  System(Body2dManifold::Instance()), 
  d(.3, .1), I(1,1,1), force(force)
{
}


double Body2d::Step(M3V3d& xb, double t, const M3V3d& xa,
                    const Vector3d& u, double h, const VectorXd* p,
                    Matrix6d *A, Matrix63d *B, Matrix<double, 6, Dynamic> *C) {
  Matrix36d fx;  // force Jacobian with respect to x
  Matrix3d fu;    // force Jacobian with respect to u
  Matrix<double, 3, Dynamic> fp;  // force Jacobian with respect to p
  
  if (p && C)
    fp.resize(3, p->size());

  if (force) {
    Vector3d f;
    force->Set(f, xa, t, u, h, p,
               A ? &fx : 0,
               B ? &fu : 0,
               C ? &fp : 0);
    xb.second = xa.second + h*f;
  } else {
    xb.second = xa.second;
  }

  Vector3d hxib = h*xb.second;
  
  SE2 &se2 = SE2::Instance();
  
  Matrix3d dg;
  
  se2.cay(dg, hxib);
  xb.first = xa.first*dg;
  
  if (A) {
    se2.cay(dg, -hxib);
    Matrix3d Adm;
    se2.Ad(Adm, dg);
    A->topLeftCorner<3,3>() = Adm;
    se2.dcay(dg, -hxib);    
    A->topRightCorner<3,3>() = h*dg;        
    if (force) {
      A->bottomRows<3>() = h*fx;
    }
    A->bottomRightCorner<3,3>() += Matrix3d::Identity();
  }
  
  if (B) {
    // A should've been called too
    assert(A);
    B->topRows<3>() = (h*h)*(dg*fu);    
    B->bottomRows<3>() = (Vector3d(h, h, h).cwiseQuotient(I)).asDiagonal()*fu;
  }

  if (C) {
    if (p) {
      assert(C->cols() == p->size() && C->rows() == 6); 
    }
    if (force) {
      assert(A);
      C->topRows<3>() = (h*h)*(dg*fp);
      C->bottomRows<3>() = (Vector3d(h, h, h).cwiseQuotient(I)).asDiagonal()*fp;
    } else {
      C->setZero();
    }
  }

  return 1;
}
