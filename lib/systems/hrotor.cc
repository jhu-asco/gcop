#include "hrotor.h"
#include "so3.h"

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

void Hrotor::StateAndControlsToFlat(VectorXd &y, const Body3dState &x,
               const Vector4d &u) {
  SO3& so3 = SO3::Instance();

  // Flat outputs are x,y,z, and yaw
  y.resize(4);
  y.head<3>() = x.p;
  y(3) = so3.yaw(x.R);
}
void Hrotor::FlatToStateAndControls(Body3dState &x, Vector4d &u,
               const std::vector<VectorXd> &y) {
  assert(y.size() >= 3);

  VectorXd y0 = y[0];
  VectorXd y1 = y[1];
  VectorXd y2 = y[2];

  x.R.setZero();
  x.p.setZero();
  x.w.setZero();
  x.v.setZero();
  u.setZero();

  Vector3d f_thrust = m*y2.head<3>() - m*Vector3d(0,0,-9.81);
  Vector3d z_rot = f_thrust/f_thrust.norm();  
  Vector3d x_yaw(cos(y0(3)), sin(y0(3)), 0);
  Vector3d y_rot = z_rot.cross(x_yaw);
  y_rot  = y_rot/y_rot.norm();
  Vector3d x_rot = y_rot.cross(z_rot);
  //Vector3d x_rot = x_yaw - x_yaw.dot(z_rot)*z_rot;
  //x_rot = x_rot/x_rot.norm();
  //Vector3d y_rot = z_rot.cross(x_rot);

  x.R.block<3,1>(0,0) = x_rot;
  x.R.block<3,1>(0,1) = y_rot;
  x.R.block<3,1>(0,2) = z_rot;

  x.p = y0.head<3>();
  //std::cout << "x: " << x.second.head<3>().transpose() << " y: " << y0.transpose() << endl;
  x.v = y1.head<3>();
  // TODO: fill in angular velocity and controls
  u(3) = f_thrust.norm();
  x.w[2] = y1(3);
}

