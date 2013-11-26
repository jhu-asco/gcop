#include "heli.h"

using namespace gcop;
using namespace Eigen;

Heli::Heli() : 
  Body3d(Vector3d(1, 1, 2), 5, Vector3d(2.32e-2, 2.32e-2, 4.41e-2)),
  rt(.3), rb(1) {
  
  Bu.setZero();
  Bu(0,0) = 1;
  Bu(1,1) = 1;
  Bu(2,2) = 1;
  Bu(5,3) = 1;
  
  fp << 0, 0, -9.81*m;

  //  double ulb[4] = {1200, 1200, 1200, 1200};
  //  double uub[4] = {7800, 7800, 7800, 7800};
}
