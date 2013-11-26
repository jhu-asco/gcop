#include "mbscspace.h"
#include "se3.h"

using namespace gcop;

MbsCspace::MbsCspace(int nb) : Manifold<MbsState>(nb + 5) {
}
  
  
void MbsCspace::Lift(VectorXd &v, 
                     const MbsState &xa,
                     const MbsState &xb) 
{
  Matrix4d gi;
  SE3::Instance().inv(gi, xa.gs[0]);
  
  Vector6d xi;
  SE3::Instance().cayinv(xi, gi*xb.gs[0]);
  
  //v.head<6>() = xi;
  v.head(6) = xi;
  //    v.tail<6 + 2*(nb - 1)>() = xb.second - xa.second;
  v.tail(xa.r.size()) = xb.r - xa.r;
}


void MbsCspace::Retract(MbsState &xb, 
                        const MbsState &xa,
                        const VectorXd &v) 
{  
  Matrix4d dg;  
  SE3::Instance().cay(dg, v.head(6));
  xb.gs[0] = xa.gs[0]*dg;
  xb.r = xa.r + v.tail(xa.r.size());
}
