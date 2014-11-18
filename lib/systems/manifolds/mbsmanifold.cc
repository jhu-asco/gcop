#include "mbsmanifold.h"  
#include "se3.h"
#include <iostream>
#include <cmath>
  

using namespace gcop;
using namespace std;

MbsManifold::MbsManifold(int nb, bool fixed) : 
  Manifold(2*(nb -1) + 12*(!fixed)), 
  cay(false) {
  lb = MbsState(nb, fixed);
  ub = MbsState(nb, fixed);
}


void MbsManifold::Lift(VectorXd &v, 
                       const MbsState &xa,
                       const MbsState &xb) {
  
  assert(!std::isnan(xa.gs[0](0,0)));

  int nb = xa.r.size() + 1;
  
  if (xa.fixed) {
    v.head(nb - 1) = xb.r - xa.r;
    v.tail(nb - 1) = xb.dr - xa.dr;    
  } else {
    Matrix4d gi;
    SE3::Instance().inv(gi, xa.gs[0]);
    
    Vector6d xi;
    if (cay)
      SE3::Instance().cayinv(xi, gi*xb.gs[0]);
    else
      SE3::Instance().log(xi, gi*xb.gs[0]);
    
    v.head(6) = xi;
    v.segment(6, nb-1) = xb.r - xa.r;
    v.segment(nb+5, 6) = xb.vs[0] - xa.vs[0];
    v.segment(nb+11, nb-1) = xb.dr - xa.dr;
  }
}


void MbsManifold::Retract(MbsState &xb, 
                          const MbsState &xa,
                          const VectorXd &v) {

  int nb = xa.r.size() + 1;

  Matrix4d dg;    
  
  if (xa.fixed) {
    xb.r = xa.r + v.head(nb-1);
    xb.dr = xa.dr + v.tail(nb-1);
  } else {
    if (cay)
      SE3::Instance().cay(dg, v.head(6));
    else
      SE3::Instance().exp(dg, v.head(6));
    
    xb.gs[0] = xa.gs[0]*dg;
    xb.r = xa.r + v.segment(6,nb-1);
    xb.vs[0] = xa.vs[0] + v.segment(nb+5, 6);
    xb.dr = xa.dr + v.segment(nb+11, nb-1);
  }
}


void MbsManifold::dtau(MatrixXd &M, const VectorXd &v)
{
  M.setIdentity();
  if (!lb.fixed) {
    Matrix6d D;
    if (cay)
      SE3::Instance().dcay(D, v.head<6>());
    else
      SE3::Instance().tln(D, -v.head<6>());
    
    M.topLeftCorner<6,6>() = D;
  }
}

void MbsManifold::Adtau(MatrixXd &M, const VectorXd &v)
{
  M.setIdentity();
  if (!lb.fixed) {
    Matrix4d g;
    if (cay)
      SE3::Instance().cay(g, v.head<6>());
    else
      SE3::Instance().exp(g, v.head<6>());
    Matrix6d A;
    SE3::Instance().Ad(A, g);
    M.topLeftCorner<6,6>() = A;
  }
}
