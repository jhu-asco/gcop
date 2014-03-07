#include "chain1.h"

using namespace gcop;
using namespace Eigen;

Chain1::Chain1() : 
  Mbs(4, 9) {

  // no gravity
  ag << 0, 0, 0.0; 
  //ag << 0, 0, -9.81; 
  
  // structural properties
  Vector3d ds1(.3, .1, .2);
  double m = 1;
	cout<<"Mass1"<<m<<endl;
  
  double l0 = .1;
  double l1 = .05;
  double l2 = .1;
  Vector3d ds2(l0, l2, l1);
  Vector3d ds3(l0, l1, l2);

  this->links[0].ds= ds1;
  this->links[1].ds= ds3;
  this->links[2].ds= ds2;
  this->links[3].ds= ds2;

  this->links[0].m= 2;
  this->links[1].m= m;
  this->links[2].m= m;
  this->links[3].m= m;

  // base body
  //  this->g0s[0].setIdentity();
  this->pis[0] = -1;

  // first body
  this->se3.rpyxyz2g(this->joints[0].gp, Vector3d(0,0,0), Vector3d(0, 0, 0.1));
  this->se3.rpyxyz2g(this->joints[0].gc, Vector3d(0,0,0), Vector3d(0, 0, 0.05));

  this->pis[1] = 0;

  this->joints[0].a.setZero();
  this->joints[0].a[1] = -1;

  // second body
  this->se3.rpyxyz2g(this->joints[1].gp, Vector3d(0,0,0), Vector3d(0.15, 0, -0.05));
  this->se3.rpyxyz2g(this->joints[1].gc, Vector3d(0,0,0), Vector3d(0.05, 0, 0));
  this->pis[2] = 0;

  this->joints[1].a.setZero();
  this->joints[1].a[1] = -1;

	// third body
  this->se3.rpyxyz2g(this->joints[2].gp, Vector3d(0,0,0), Vector3d(-0.15, 0, -0.05));
  this->se3.rpyxyz2g(this->joints[2].gc, Vector3d(0,0,0), Vector3d(-0.05, 0, 0));
  this->pis[3] = 0;

  this->joints[2].a.setZero();
  this->joints[2].a[1] = -1;

	this->links[0].I(0) =0.01833 ;
	this->links[0].I(1) = 0.045;
	this->links[0].I(2) = 0.04;
	this->links[0].I(3) = 2;
	this->links[0].I(4) = 2;
	this->links[0].I(5) = 2;
  Body3d<>::Compute(this->links[1].I, m, ds3);
  Body3d<>::Compute(this->links[2].I, m, ds2);
  Body3d<>::Compute(this->links[3].I, m, ds2);

	this->debug = true;
  
  this->Init();
}

void Chain1::Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
                  MatrixXd *A, MatrixXd *B) 
{
  f = u;
  if (A)
    A->setZero();
  if (B)
    B->setIdentity();
}
