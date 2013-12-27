#include "chain.h"

using namespace gcop;
using namespace Eigen;

Chain::Chain() : 
  Mbs(3, 8) {

  // no gravity
  //ag << 0, 0, -9.81; 
  
  // structural properties
  Vector3d ds(.3, .1, .1);
  double m = 10;
	cout<<"Mass1"<<m<<endl;
  
  double l0 = .3;
  double l1 = .3;
  double l2 = .3;

  this->links[0].ds= ds;
  this->links[1].ds= ds;
  this->links[2].ds= ds;

  this->links[0].m= m;
  this->links[1].m= m;
  this->links[2].m= m;

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

void Chain::Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
                  MatrixXd *A, MatrixXd *B) 
{
  f = u;
  if (A)
    A->setZero();
  if (B)
    B->setIdentity();
}
