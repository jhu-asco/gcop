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

void Hrotor::StateAndControlsToFlatAndDerivatives(vector<VectorXd> &y, 
                                          const Body3dState &x, const std::vector<Vector4d> &u)
{
  assert(u.size() >= 2);
  SO3& so3 = SO3::Instance();
  y.clear();
  y.resize(5);
  for(int i = 0; i < y.size(); i++)
  {
    y.at(i).resize(4);
  }

  // 0th 
  y.at(0).head<3>() = x.p; 
  y.at(0)(3) = so3.yaw(x.R);

  // 1st
  y.at(1).head<3>() = x.v;
  y.at(1)(3) = x.w(2)/x.R(2,2);

  // 2nd
  y.at(2).head<3>() = (x.R*Eigen::Vector3d(0,0,u.at(0)(3)) + m*Vector3d(0,0,-9.81))/m;
  y.at(2)(3) = 0;

  // 3rd
  Eigen::Matrix3d what;
  so3.hat(what, x.w);
  y.at(3).head<3>() = (x.R*what*Eigen::Vector3d(0,0,u.at(0)(3)) 
                      + x.R*Eigen::Vector3d(0,0,u.at(1)(3)))/m;
  y.at(3)(3) = 0;

  // 4th
  y.at(4).setZero();
}

void Hrotor::StateAndControlsToFlat(VectorXd &y, const Body3dState &x,
               const Vector4d &u) {
  SO3& so3 = SO3::Instance();

  // Flat outputs are x,y,z, and yaw
  y.resize(4);
  y.head<3>() = x.p;
  y(3) = so3.yaw(x.R);
}
void Hrotor::FlatToStateAndControls(Body3dState &x, std::vector<Vector4d> &u,
               const std::vector<VectorXd> &y) {
  assert(y.size() >= 5);
  SO3& so3 = SO3::Instance();
  u.resize(2);
  //SO3& so3 = SO3::Instance();

  VectorXd y0 = y[0];
  VectorXd y1 = y[1];
  VectorXd y2 = y[2];
  VectorXd y3 = y[3];
  VectorXd y4 = y[4];

  x.Clear();
  u.at(0).setZero();
  u.at(1).setZero();

  Vector3d f_thrust = m*y2.head<3>() - m*Vector3d(0,0,-9.81);
  Vector3d z_rot = f_thrust/f_thrust.norm();
  double x3 = (-z_rot(0)*cos(y0(3)) - z_rot(1)*sin(y0(3)))/z_rot(2);  
  Vector3d x_yaw(cos(y0(3)), sin(y0(3)), x3);
  Vector3d x_rot = x_yaw.normalized();
  Vector3d y_rot = z_rot.cross(x_rot);

  //Vector3d y_rot = z_rot.cross(x_yaw);
  //y_rot  = y_rot/y_rot.norm();
  //Vector3d x_rot = y_rot.cross(z_rot);
  //Vector3d x_rot = x_yaw - x_yaw.dot(z_rot)*z_rot;
  //x_rot = x_rot/x_rot.norm();
  //Vector3d y_rot = z_rot.cross(x_rot);

  x.R.block<3,1>(0,0) = x_rot;
  x.R.block<3,1>(0,1) = y_rot;
  x.R.block<3,1>(0,2) = z_rot;
  //std::cout << x.R.determinant() << " " << x.R.transpose()*x.R << std::endl;

  x.p = y0.head<3>();
  //std::cout << "x: " << x.second.head<3>().transpose() << " y: " << y0.transpose() << endl;
  x.v = y1.head<3>();
  u.at(0)(3) = f_thrust.norm();
  x.w[0] = -m*y_rot.transpose().dot(y3.head<3>())/u.at(0)(3);
  x.w[1] = m*x_rot.transpose().dot(y3.head<3>())/u.at(0)(3);
  x.w[2] = y1(3)*z_rot(2);
  
  u.at(1)(3) = m*z_rot.dot(y3.head<3>());

  // torque controls
  Eigen::Matrix3d what;
  so3.hat(what, x.w);
  
  Vector3d wdot;
  wdot(2) = y2(3)*z_rot(2) + y1(3)*(x.R*what)(2,2);
  wdot(1) = (m*x_rot.dot(y4.head<3>()) - x.w(0)*x.w(2)* u.at(0)(3) 
            - 2*x.w(1)*u.at(1)(3))/u.at(0)(3);
  wdot(0) = -(m*y_rot.dot(y4.head<3>()) - x.w(1)*x.w(2)* u.at(0)(3) 
            + 2*x.w(0)*u.at(1)(3))/u.at(0)(3);
  u.at(0).head<3>() = J.asDiagonal()*wdot-(J.asDiagonal()*x.w).cross(x.w);



  /*
  std::vector<VectorXd> ycheck;
  StateAndControlsToFlatAndDerivatives(ycheck, x, u);
  for(int i = 0; i < ycheck.size(); i++)
  {
    std::cout << "der " << i << "= " << ycheck[i].transpose() << " <> " << y[i].transpose() << std::endl;
  }
  */
}

