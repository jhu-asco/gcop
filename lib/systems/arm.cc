#include "arm.h"
#include <stdio.h>
#include <cmath>

using namespace gcop;

Arm::Arm() : l1(.2), l2(.15), x1(.03)
{
}


Arm::~Arm()
{
}

/*
  FK corresponds to the end-eff transformation

T = se3_exp([0; a1; 0; 0;0;0])*...
    se3_exp([0; 0; 0; x1; 0; l1])*...
    se3_exp([0; a2; 0; 0; 0; 0])*...
    se3_exp([0; 0; 0; 0; 0; l2]);
*/
bool Arm::Fk(double p[3], const double a[3])
{
  const double &a0 = a[0]; // yaw
  const double &a1 = a[1]; // first joint
  const double &a2 = a[2]; // second joint
  
  double xr = l2*(cos(a1)*sin(a2) + cos(a2)*sin(a1)) + x1*cos(a1) + l1*sin(a1);
  double zr = l1*cos(a1) - l2*(sin(a1)*sin(a2) - cos(a1)*cos(a2)) - x1*sin(a1);

  p[0] = cos(a0)*xr;
  p[1] = sin(a0)*xr;
  p[2] = zr;

  return true;
}


bool Arm::Ik(double a[2][3], const double p[3])
{
  // it easier to derive IK in terms of polar coordinates
  double r = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance to end-eff

  double yaw = atan2(p[1], p[0]);
  double xr = cos(yaw)*p[0] + sin(yaw)*p[1];
  double yr = -sin(yaw)*p[0] + cos(yaw)*p[1];

  double c = atan2(p[2], xr);                         // angle to end-eff

  a[0][0] = yaw;     // yaw (assume vehicle is pointing in the direction as arm is pointing)
  a[1][0] = yaw;

  double l1_2 = l1*l1;
  double l1_4 = l1_2*l1_2;
  double l2_2 = l2*l2;
  double l2_4 = l2_2*l2_2;
  double r_2 = r*r;
  double r_4 = r_2*r_2;
  double x1_2 = x1*x1;
  double x1_4 = x1_2*x1_2;

  double C = - l1_4 + 2*l1_2*l2_2 + 2*l1_2*r_2 - 2*l1_2*x1_2 - l2_4 + 2*l2_2*r_2 + 2*l2_2*x1_2 - r_4 + 2*r_2*x1_2 - x1_4;
  // outside of manipulator workspace
  if (C < 0)
    return false;

  double ks[2][2];

  // Grobner basis
  ks[0][0] = -(sqrt(C) - 2*l1*r*cos(c) + 2*r*x1*sin(c))/(l1_2 + 2*sin(c)*l1*r - l2_2 + r_2 + 2*cos(c)*r*x1 + x1_2);
  ks[1][0] = (sqrt(C) + 2*l1*r*cos(c) - 2*r*x1*sin(c))/(l1_2 + 2*sin(c)*l1*r - l2_2 + r_2 + 2*cos(c)*r*x1 + x1_2);

  ks[0][1] = -(2*l2*x1 + sqrt(C))/(l1_2 - 2*l1*l2 + l2_2 - r_2 + x1_2);
  ks[1][1] = -(2*l2*x1 - sqrt(C))/(l1_2 - 2*l1*l2 + l2_2 - r_2 + x1_2);
  
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) { 
      double &k = ks[i][j];
      a[i][1 + j] = atan2(2*k/(1+k*k), (1-k*k)/(1+k*k));
    }
  }

  return true;
}
