#include "chain.h"

using namespace gcop;
using namespace Eigen;

Chain::Chain(int nb, bool fixed) : 
  Mbs(nb, nb - 1 + 6*(!fixed), fixed) {

  if (fixed)
    basetype = FIXEDBASE;
  
  // no gravity
  //ag << 0, 0, -9.81; 
  
  // structural properties
  Vector3d ds(.3, .1, .1);
  double m = 1;  
  double l0 = .3;
  double l1 = .3;
  double l2 = .3;

  this->pis[0] = -1;

  for (int i = 0; i < nb; ++i) {
    this->links[i].ds = ds;
    this->links[i].m = m;
    Body3d<>::Compute(this->links[i].I, m, ds);    
    this->pis[i] = i - 1;
  }    

  for (int i = 0; i < nb - 1; ++i) {
    this->se3.rpyxyz2g(this->joints[i].gp, Vector3d(0,0,0), Vector3d(l0/2, 0, 0));
    this->se3.rpyxyz2g(this->joints[i].gc, Vector3d(0,0,0), Vector3d(-l1/2, 0, 0));
    
    this->joints[i].a.setZero();
    this->joints[i].a[2] = 1;
    this->joints[i].lower = -3;
    this->joints[i].upper = 3;    
    X.lb.r[i] = this->joints[i].lower;
    X.ub.r[i] = this->joints[i].upper;      

    X.lb.dr[i] = -10;
    X.ub.dr[i] = 10;
  }

  U.bnd = true;
  if (!fixed) {
    for (int i = 0; i < 3; ++i) {
      X.lb.gs[0](i,3) = -10;
      X.ub.gs[0](i,3) = 10;
    }
    X.lb.vs[0] << -10, -10, -10, -100, -100, -100;
    X.ub.vs[0] << 10, 10, 10, 100, 100, 100;    

    for (int i = 0; i < 3; ++i) {
      U.lb[i] = -2;
      U.ub[i] = 2;
    }
    for (int i = 3; i < 6; ++i) {
      U.lb[i] = -10;
      U.ub[i] = 10;
    }
    for (int i = 6; i < nb+5; ++i) {
      U.lb[i] = -1;
      U.ub[i] = 1;
    }
  } else {
    for (int i = 0; i < nb-1; ++i) {
      U.lb[i] = -1;
      U.ub[i] = 1;
    }
  }

  this->Init();
}

/*
void Chain::Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
                  MatrixXd *A, MatrixXd *B) 
{
  if (A)
    A->setZero();
  if (B)
    B->setIdentity();
}
*/
