#include <limits>
#include "airplanecost.h"
#include "gunicyclemanifold.h"
#include "se2.h"
#include <iostream>
#include <assert.h>

using namespace gcop;
using namespace Eigen;

AirplaneCost::AirplaneCost(double tf, const M3V2d &xf, int N) : 
  LqCost(GunicycleManifold::Instance(), tf, xf), gds(N+1, Normal(2))
{
  Q(0,0) = .01;
  Q(1,1) = .01;
  Q(2,2) = .01;
  Q(3,3) = .5;
  Q(4,4) = .5;
  
  Qf(0,0) = 5;
  Qf(1,1) = 5;
  Qf(2,2) = 5;
  Qf(3,3) = 10;
  Qf(4,4) = 10;
  
  R(0,0) = 1;
  R(1,1) = 2;
  
  // stdev
  double s = 1;
  for (int i =0; i <= N; ++i) {
    gds[i].P << 1/(s*s), 0, 0, 1/(s*s);
    gds[i].mu << 5, ((N/2-i)*10.0/N);
    gds[i].Update();
  }
  n = 10;
}


void AirplaneCost::SetObstacle(const double s, double n)
{
  for (int i =0; i < gds.size(); ++i) {
    gds[i].P << s*s, 0, 0, s*s;
    //    gd.mu = po;
    gds[i].Update();
  }
  this->n = n;    
}


double AirplaneCost::L(double t, const M3V2d& x, const Vector2d& u, double h,
                       Vector5d *Lx, Matrix5d* Lxx,
                       Vector2d *Lu, Matrix2d* Luu,
                       Matrix52d *Lxu) 
{
  double c0 = LqCost::L(t, x, u, h, Lx, Lxx, Lu, Luu, Lxu);
    
  if (t > tf - 1e-10) {
    return c0;

  } else {
    int i = t/h;
    assert(i <= gds.size());
    
    Vector2d p = x.first.topRightCorner<2,1>();
    Matrix2d R = x.first.topLeftCorner<2,2>();
    
    Vector2d dp = p - gds[i].mu;
    double c = n*gds[i].L(p);
    
    if (Lx) {
      (*Lx).segment<2>(1) += -h*c*R.transpose()*gds[i].Pinv*dp;
    }
    
    if (Lxx) {
      (*Lxx).block<2,2>(1,1) += h*c*R.transpose()*gds[i].Pinv*(dp*dp.transpose()*gds[i].Pinv - Matrix2d::Identity())*R;      
    }
    
    return c0 + h*c;
  }
}
