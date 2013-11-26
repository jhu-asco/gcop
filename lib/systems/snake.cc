#include "snake.h"

using namespace gcop;
using namespace Eigen;

Snake::Snake() : 
  Mbs(3, 6) {

  // structural properties
  Vector3d ds(.3, .1, .1);
  double m = 1;
  
  double l0 = .3;
  double l1 = .3;
  double l2 = .3;

  this->links[0].ds= ds;
  this->links[1].ds= ds;
  this->links[2].ds= ds;

  // base body
  //  this->g0s[0].setIdentity();
  this->pis[0] = -1;

  // first body
  this->se3.rpyxyz2g(this->joints[0].gp, Vector3d(0,0,0), Vector3d(l0/2, 0, 0));
  this->se3.rpyxyz2g(this->joints[0].gc, Vector3d(0,0,0), Vector3d(-l1/2, 0, 0));

  this->pis[1] = 0;

  this->joints[0].a.setZero();
  this->joints[0].a[2] = 1;

  // second body
  this->se3.rpyxyz2g(this->joints[1].gp, Vector3d(0,0,0), Vector3d(l1/2, 0, 0));
  this->se3.rpyxyz2g(this->joints[1].gc, Vector3d(0,0,0), Vector3d(-l2/2, 0, 0));
  this->pis[2] = 1;

  this->joints[1].a.setZero();
  this->joints[1].a[2] = 1;

  Body3d<>::Compute(this->links[0].I, m, ds);
  Body3d<>::Compute(this->links[1].I, m, ds);
  Body3d<>::Compute(this->links[2].I, m, ds);  
  
  this->Init();
}

void Snake::Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
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
  f[6] = u[4];
  f[7] = u[5];

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
