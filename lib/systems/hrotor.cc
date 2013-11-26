#include "hrotor.h"

using namespace gcop;
using namespace Eigen;

Hrotor::Hrotor() : 
  Body3d(Vector3d(.5, .5, .1), .5, Vector3d(2.32e-3, 2.32e-3, 4.41e-3)),
  l(0.171), r(0.099), 
  kt(6.107e-8), km(2.7207e-9), mm(20) {
  
  Bu.setZero();
  Bu(0,0) = 1;
  Bu(1,1) = 1;
  Bu(2,2) = 1;
  Bu(5,3) = 1;
  
  fp << 0, 0, -9.81*m;

  //  double ulb[4] = {1200, 1200, 1200, 1200};
  //  double uub[4] = {7800, 7800, 7800, 7800};
}

