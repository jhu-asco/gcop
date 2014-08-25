#include "airbot.h"

using namespace gcop;
using namespace Eigen;

Airbot::Airbot() : 
  Mbs(7, 10) {

  basetype = AIRBASE;

//  ag << 0, 0, 0.0; 
  hrotor.l = 1.5*hrotor.l;
  hrotor.r = 1.5*hrotor.r;

  double z0 = .1;
  double l1 = .1;
  double l2 = .3;
  double l3 = .3;

  // structural properties
  double m0 = 1;
  double m1 = .2;
  double m2 = .3;
  double m3 = .2;

  bool limits = true;
  
  this->links[1].ds[0] = .05; this->links[1].ds[1] = l1; this->links[1].ds[2] = .05;
  this->links[2].ds[0] = .05; this->links[2].ds[1] = .05; this->links[2].ds[2] = l2;
  this->links[3].ds[0] = .05; this->links[3].ds[1] = .05; this->links[3].ds[2] = l3;

  this->links[4].ds[0] = .05; this->links[4].ds[1] = l1; this->links[4].ds[2] = .05;
  this->links[5].ds[0] = .05; this->links[5].ds[1] = .05; this->links[5].ds[2] = l2;
  this->links[6].ds[0] = .05; this->links[6].ds[1] = .05; this->links[6].ds[2] = l3;

  
  // base body
  this->pis[0] = -1;

  // body #1
  this->se3.rpyxyz2g(this->joints[0].gp, Vector3d(0,0,0), Vector3d(0, l1, -z0));
  this->se3.rpyxyz2g(this->joints[0].gc, Vector3d(0,0,0), Vector3d(0, 0, 0));

  this->pis[1] = 0;
  
  this->joints[0].a.setZero();
  this->joints[0].a[2] = 1;
  if (limits) {
    this->joints[0].lower = -3;
    this->joints[0].upper = 3;
  }
  // body #2
  this->se3.rpyxyz2g(this->joints[1].gp, Vector3d(0,0,0), Vector3d(0, l1/2,0));
  this->se3.rpyxyz2g(this->joints[1].gc, Vector3d(0,0,0), Vector3d(0,0, l2/2));

  this->pis[2] = 1;
  
  this->joints[1].a.setZero();
  this->joints[1].a[1] = 1;
  if (limits) {
    this->joints[1].lower = -3;
    this->joints[1].upper = 3;
  }

  // body #3
  this->se3.rpyxyz2g(this->joints[2].gp, Vector3d(0,0,0), Vector3d(0,0,-l2/2));
  this->se3.rpyxyz2g(this->joints[2].gc, Vector3d(0,0,0), Vector3d(0,0, l3/2));

  this->pis[3] = 2;

  this->joints[2].a.setZero();
  this->joints[2].a[1] = 1;
  if (limits) {
    this->joints[2].lower = -3;
    this->joints[2].upper = 3;
  }

  // body #4
  this->se3.rpyxyz2g(this->joints[3].gp, Vector3d(0,0,0), Vector3d(0, -l1, -z0));
  this->se3.rpyxyz2g(this->joints[3].gc, Vector3d(0,0,0), Vector3d(0, 0, 0));

  this->pis[4] = 0;
  
  this->joints[3].a.setZero();
  this->joints[3].a[2] = 1;
  if (limits) {
    this->joints[3].lower = -3;
    this->joints[3].upper = 3;
  }

  // body #5
  this->se3.rpyxyz2g(this->joints[4].gp, Vector3d(0,0,0), Vector3d(0, -l1/2,0));
  this->se3.rpyxyz2g(this->joints[4].gc, Vector3d(0,0,0), Vector3d(0, 0, l2/2));

  this->pis[5] = 4;
  
  this->joints[4].a.setZero();
  this->joints[4].a[1] = 1;
  if (limits) {
    this->joints[4].lower = -3;
    this->joints[4].upper = 3;
  }

  // body #6
  this->se3.rpyxyz2g(this->joints[5].gp, Vector3d(0,0,0), Vector3d(0,0,-l2/2));
  this->se3.rpyxyz2g(this->joints[5].gc, Vector3d(0,0,0), Vector3d(0,0, l3/2));

  this->pis[6] = 5;

  this->joints[5].a.setZero();
  this->joints[5].a[1] = 1;
  if (limits) {
    this->joints[5].lower = -3;
    this->joints[5].upper = 3;
  }

  links[0].I(0) = m0*4.32e-3;
  links[0].I(1) = m0*4.32e-3;
  links[0].I(2) = m0*8.41e-3;
  links[0].I(3) = m0;
  links[0].I(4) = m0;
  links[0].I(5) = m0;

  Body3d<>::Compute(this->links[1].I, m1, this->links[1].ds);
  Body3d<>::Compute(this->links[2].I, m2, this->links[2].ds);  
  Body3d<>::Compute(this->links[3].I, m3, this->links[3].ds);
  Body3d<>::Compute(this->links[4].I, m1, this->links[4].ds);  
  Body3d<>::Compute(this->links[5].I, m2, this->links[5].ds);
  Body3d<>::Compute(this->links[6].I, m3, this->links[6].ds);  

  this->links[0].m = m0;
  this->links[1].m = m1;
  this->links[2].m = m2;
  this->links[3].m = m3;
  this->links[4].m = m1;
  this->links[5].m = m2;
  this->links[6].m = m3;

  for (int i = 0; i < nb - 1; ++i) {
    X.lb.r[i] = this->joints[i].lower;
    X.ub.r[i] = this->joints[i].upper;      
    X.lb.dr[i] = -10;
    X.ub.dr[i] = 10;
  }

  U.bnd = true;

  for (int i = 0; i < 3; ++i) {
    X.lb.gs[0](i,3) = -10;
    X.ub.gs[0](i,3) = 10;
  }

  X.lb.vs[0] << -10, -10, -10, -100, -100, -100;
  X.ub.vs[0] << 10, 10, 10, 100, 100, 100;    
  
  // max torques
  for (int i = 0; i < 3; ++i) {
    U.lb[i] = -1;
    U.ub[i] = 1;
  }
  U.lb[3] = 0;
  U.ub[3] = 50; // max lift

  for (int i = 4; i < U.n; ++i) {
    U.lb[i] = -1;
    U.ub[i] = 1;
  }
  
  this->Init();
}

void Airbot::Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
                   MatrixXd *A, MatrixXd *B) 
{
  assert(f.size() == 12);
  assert(u.size() == 10);

  f[0] = u[0];
  f[1] = u[1];
  f[2] = u[2];
  f[3] = 0;
  f[4] = 0;
  f[5] = u[3];
  
  f.tail(nb-1) = u.tail(nb-1) - damping.cwiseProduct(x.dr);
  
  //  if (A) 
  //    A->setZero();

  if (B) {
    B->setZero();
    (*B)(0,0) = 1;
    (*B)(1,1) = 1;
    (*B)(2,2) = 1;
    (*B)(5,3) = 1;
    (*B)(6,4) = 1;
    (*B)(7,5) = 1;
    (*B)(8,6) = 1;
    (*B)(9,7) = 1;
    (*B)(10,8) = 1;
    (*B)(11,9) = 1;
  }
}
