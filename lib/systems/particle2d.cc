#include "particle2d.h"
#include <limits>
#include "rn.h"

using namespace gcop;
using namespace Eigen;
 
Particle2d::Particle2d() : System(Rn<4>::Instance(), Rn<2>::Instance()), m(1), r(.1)
{
}


double Particle2d::Step(Vector4d& xb, double t, const Vector4d& xa,
                        const Vector2d& u, double h,
                        Matrix4d *A, Matrix<double, 4, 2> *B) {

  xb.head(2) = xa.head(2) + h*xa.tail(2);
  xb.tail(2) = xa.tail(2) + (h/m)*u;
  
  if (A) {
    (*A).block<2,2>(0,0).setIdentity();
    (*A)(0,2) = h, (*A)(0, 3) =  0;
    (*A)(1,2) = 0, (*A)(1, 3) =  h;
    (*A).block<2,2>(2,0).setZero();
    (*A).block<2,2>(2,2).setIdentity();
  }
  
  if (B) {
    (*B).block<2,2>(0,0).setZero();
    (*B)(2,0) = h/m, (*B)(2, 1) =  0;
    (*B)(3,0) = 0, (*B)(3, 1) =  h/m;      
  }  
  return 1;
}
