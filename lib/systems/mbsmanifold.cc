#include "mbsmanifold.h"  
#include "se3.h"
#include <iostream>
#include <cmath>
  

using namespace gcop;
using namespace std;

MbsManifold::MbsManifold(int nb) : Manifold(MBS_DIM(nb)) {
}


void MbsManifold::Lift(VectorXd &v, 
                       const MbsState &xa,
                       const MbsState &xb) {
  
  assert(!std::isnan(xa.gs[0](0,0)));

  int nb = xa.r.size() + 1;

  Matrix4d gi;
  SE3::Instance().inv(gi, xa.gs[0]);
  
  Vector6d xi;
  SE3::Instance().cayinv(xi, gi*xb.gs[0]);

  /*
  cout << "xa.gs[0]=" << xa.gs[0] << endl;
  cout << "xa.vs[0]=" << xa.vs[0] << endl;
  cout << "xa.r=" << xa.r << endl;
  cout << "xa.dr=" << xa.dr << endl;

  cout << "xb.gs[0]=" << xb.gs[0] << endl;
  cout << "xb.vs[0]=" << xb.vs[0] << endl;
  cout << "xb.r=" << xb.r << endl;
  cout << "xb.dr=" << xb.dr << endl;
  */
  
  //v.head<6>() = xi;
  v.head(6) = xi;
  //    v.tail<6 + 2*(nb - 1)>() = xb.second - xa.second;
  v.segment(6, nb-1) = xb.r - xa.r;
  v.segment(nb+5, 6) = xb.vs[0] - xa.vs[0];    
  v.segment(nb+11, nb-1) = xb.dr - xa.dr;
}


void MbsManifold::Retract(MbsState &xb, 
                          const MbsState &xa,
                          const VectorXd &v) {

  int nb = xa.r.size() + 1;

  Matrix4d dg;    

  SE3::Instance().cay(dg, v.head(6));
  xb.gs[0] = xa.gs[0]*dg;
  xb.r = xa.r + v.segment(6,nb-1);
  xb.vs[0] = xa.vs[0] + v.segment(nb+5, 6);
  xb.dr = xa.dr + v.segment(nb+11, nb-1);
}


void MbsManifold::dtau(MatrixXd &M, const VectorXd &v)
{
  M.setIdentity();
  Matrix6d D;
  SE3::Instance().dcay(D, v.head<6>());
  M.topLeftCorner<6,6>() = D;
}

void MbsManifold::Adtau(MatrixXd &M, const VectorXd &v)
{
  M.setIdentity();
  Matrix4d g;
  SE3::Instance().cay(g, v.head<6>());
  Matrix6d A;
  SE3::Instance().Ad(A, g);
  M.topLeftCorner<6,6>() = A;
}
