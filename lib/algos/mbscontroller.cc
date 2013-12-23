#include "mbscontroller.h"

using namespace gcop;
using namespace std;
using namespace Eigen;

MbsController::MbsController(const Mbs &sys, 
                             MbsState *xd, 
                             VectorXd *ad) : 
  Controller(),
  sys(sys), xd(xd), ad(ad), Kp(5 + sys.nb), Kd(5 + sys.nb)
{
  Kp.setOnes();
  Kd.setOnes();
}

MbsController::~MbsController()
{

}

void MbsController::Set(VectorXd &u, double t, const MbsState &x)
{
  // error in base body configuration
  Vector6d ge;

  // full error in configuration
  VectorXd e(5 + sys.nb); 

  if (xd) {
    Matrix4d gi;
    SE3::Instance().inv(gi, xd->gs[0]);
    SE3::Instance().log(ge, gi*x.gs[0]);
    
    e.head(6) = ge;
    e.tail(sys.nb-1) = x.r - xd->r;    
    
  } else {
    SE3::Instance().log(ge, x.gs[0]);
    e.head(6) = ge;
    e.tail(sys.nb-1) = x.r;
  }
  
  // full error in velocity
  VectorXd de(5 + sys.nb); 

  if (xd) {
    de.head(6) = x.vs[0] - xd->vs[0];
    de.tail(sys.nb - 1) = x.dr - xd->dr;
  } else {
    de.head(6) = x.vs[0];
    de.tail(sys.nb - 1) = x.dr;
  }

  // desired change in velocity
  VectorXd dv = -Kp.cwiseProduct(e) - Kd.cwiseProduct(de);
  // add desired acceleration if provided
  if (ad)
    dv = dv + *ad;

  MatrixXd M(sys.nb + 5, sys.nb + 5);
  sys.Mass(M, x);

  VectorXd b(sys.nb + 5);  // bias forces
  
  // compute bias b (bias means that M*dv + b = u)
  sys.Bias(b, t, x);

  // assume there is no control matrix transformation for now
  u = M*dv + b;
}

