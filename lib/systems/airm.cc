#include "airm.h"

using namespace gcop;
using namespace Eigen;

Airm::Airm() : 
  Mbs(3,6) {

  double l0 = .0;
  double l1 = .3;
  double l2 = .3;

  // structural properties
  double m0 = 1;
  double m1 = .2;
  double m2 = .2;
  
  this->links[1].ds[0] = .05; this->links[1].ds[1] = .05; this->links[1].ds[2] = l1;
  this->links[2].ds[0] = .05; this->links[2].ds[1] = .05; this->links[2].ds[2] = l2;

  // base body
  this->pis[0] = -1;

  // first body
  this->se3.rpyxyz2g(this->joints[0].gp, Vector3d(0,0,0), Vector3d(0,0,0));
  this->se3.rpyxyz2g(this->joints[0].gc, Vector3d(0,0,0), Vector3d(0,0, l1/2));

  this->pis[1] = 0;
  
  this->joints[0].a.setZero();
  this->joints[0].a[1] = 1;

  // first body
  this->se3.rpyxyz2g(this->joints[1].gp, Vector3d(0,0,0), Vector3d(0,0,-l1/2));
  this->se3.rpyxyz2g(this->joints[1].gc, Vector3d(0,0,0), Vector3d(0,0, l2/2));

  this->pis[2] = 1;

  this->joints[1].a.setZero();
  this->joints[1].a[1] = 1;

  links[0].I(0) = m0*4.32e-3;
  links[0].I(1) = m0*4.32e-3;
  links[0].I(2) = m0*8.41e-3;
  links[0].I(3) = m0;
  links[0].I(4) = m0;
  links[0].I(5) = m0;

  Body3d<>::Compute(this->links[1].I, m1, this->links[1].ds);
  Body3d<>::Compute(this->links[2].I, m2, this->links[1].ds);  
  
  this->Init();
}

void Airm::Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
                 MatrixXd *A, MatrixXd *B) 
{
  assert(f.size() == 8);
  assert(u.size() == 6);

  f[0] = u[0];
  f[1] = u[1];
  f[2] = u[2];
  f[3] = 0;
  f[4] = 0;
  f[5] = u[3];
  f[6] = u[4] - .01*x.dr[0];
  f[7] = u[5] - .01*x.dr[1];

  if (A)
    A->setZero();

  if (B) {
    B->setZero();
    (*B)(0,0) = 1;
    (*B)(1,1) = 1;
    (*B)(2,2) = 1;
    (*B)(5,3) = 1;
    (*B)(6,4) = 1;
    (*B)(7,5) = 1;
  }
}
