#include "gavoidcontroller.h"
#include "utils.h"

using namespace gcop;
using namespace std;
using namespace Eigen;

GavoidController::GavoidController(Body3dConstraint &con) :
  con(con)
{
  sr = 10;
  k = 1;
  kb = 0;
}

GavoidController::~GavoidController()
{

}

void GavoidController::Set(Vector3d &u, double t, const Body3dState &x)
{
  Matrix<double, 1, 12> dgdx;
  Matrix<double, 1,1> d; // distance to obstacle
  con(d, t, x, u, 0, &dgdx);
  
  Vector3d N = dgdx.segment<3>(3); // unit vector to obstacle

  Vector3d S(0,0,0);  // axis of rotation
  
  double od = -d[0];   // in GCOP constraints are defined as c(x,u) <= 0 so need to invert
  //  cout << "od=" << od << " N=" << N << endl;
      
  //  Vector3d &p = x.second.head<3>(); // position
  Vector3d v = x.second.tail<3>(); // velocity

  double tol = 1e-12;
  
  Vector3d fb(0,0,0);  // optional breaking force

  if (od < this->sr) {
    double vn = v.norm();
    if (vn < tol)
      return;
    
    Vector3d Sv = N.cross(v/vn);
    
    double b = 0;
    if (Sv.norm() < tol) {
      Vector3d rv(3);
      rv[0] = RND; 
      rv[1] = RND;
      rv[2] = RND;
      rv = rv/rv.norm();
      Sv = rv.cross(N);           
    } else {
      double a = asin(Sv.norm());
      b = SIGN(a)*acos(N.dot(v/vn));
    }
    //  cout << "b=" << b << endl;
    if (fabs(b) < M_PI/2) {
      S = S + this->k/od*Sv/Sv.norm();
    }
    fb = fb + (this->kb/od)*v;
  }

  // cout << "S=" << S << " fS=" << cross(S, sys.stvel ? v : R*v) << endl;
  u = S.cross(v) - fb;  
}
