#ifndef GCOP_KINRCCARPATH_H
#define GCOP_KINRCCARPATH_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include "se3.h"
#include <iostream>
#include "utils.h"
#include "kinrccar.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  
  /**
   * Kinematic RC Car Pose-track in 3D. Tracks odometry, controls, and landmarks.
   */
  class KinRccarPath {
  public:
    /**
     * Pose graph (using given features)
     * @param sys body3d system
     * @param odometry is odometry available
     * @param forces use uncertain forces in the cost function
     */
    KinRccarPath(KinRccar &sys, bool wheel_odometry = false, bool yaw_odometry = false, bool forces = false, bool estimate_length = false);

    virtual void AddControl(const Vector2d &u, const Vector2d& v, double h);
    virtual void AddObservation(const vector< pair<int, Vector3d> > &z, Matrix4d cam_transform = MatrixXd::Identity(4,4));

    KinRccar &sys;     ///< system

    bool wheel_odometry;
    bool yaw_odometry;
    bool extforce;         ///< is there a constant external force that need to be estimated
    bool forces;           ///< is there a constant external force that need to be estimated
    bool estimate_length;  ///< should the length of the car be tracked as a parameter
    
    vector<double> ts;     ///< sequence of times (N+1 vector)

    vector<Matrix4d> xs;      ///< sequence of states x=(g,v) (N+1 vector)
    vector<Vector2d> us;   ///< sequence of inputs (N vector)
    vector<Vector2d> vs;   ///< sequence of odometry (N vector)
    
    vector<Matrix4d> xos;     ///< unoptimized trajectory
    vector<Vector2d> uos;     ///< unoptimized trajectory

    VectorXd p;            ///< observed landmark positions (2*m vector) and external force parameters (2 vector)    

    vector< vector<pair<int, Vector3d> > > zs;  ///< vector of vectors of landmark observations

    double cp;             ///< noise covariance of poses (assume spherical)
    Vector2d cv;           ///< noise covariance of odometry/gyro (assume diagonal)
    Vector2d cw;           ///< noise covariance of control forces (assume diagonal)

  };

KinRccarPath::KinRccarPath(KinRccar &sys, bool wheel_odometry, bool yaw_odometry, bool forces, bool estimate_length) : 
  sys(sys), wheel_odometry(wheel_odometry), yaw_odometry(yaw_odometry), forces(forces), extforce(false), estimate_length(estimate_length),  p(estimate_length + extforce*3)
{
  Matrix4d x = MatrixXd::Identity(4,4);
  xs.push_back(x);
  xos.push_back(x);
  ts.push_back(0);

  cw = Vector2d(.01, pow(1.5*M_PI/180.,2));
  cv = Vector2d(.01, pow(2.5*M_PI/180.,2));
  cp = .01;

  if(estimate_length)
  {
    p(0) = .25;    
  }
}

void KinRccarPath::AddControl(const Vector2d &u, const Vector2d& v, double h)
{
  Eigen::Matrix4d xn;
  double t = ts.back();

  sys.Step(xn, t, xs.back(), u, h);
  xs.push_back(xn);
  us.push_back(u);

  sys.Step(xn, t, xos.back(), u, h);
  xos.push_back(xn);
  uos.push_back(u);

  vs.push_back(v);
  ts.push_back(t + h);
  zs.push_back(vector<pair<int, Vector3d> >());
}

void KinRccarPath::AddObservation(const vector< pair<int, Vector3d> > &z, Matrix4d cam_transform)
{
  for(int i = 0; i < z.size(); i++)
  {
    Vector4d z_h(z[i].second(0), z[i].second(1), z[i].second(2), 1);
    z_h = cam_transform*z_h;
    pair<int, Vector3d> zi(z[i].first, Vector3d(z_h(0), z_h(1), z_h(2)));
    if(3*z[i].first >= p.size()-estimate_length-extforce*3)
    {
        p.conservativeResize(p.size() + 3);
        const Matrix3d &Rn = xs.back().block<3,3>(0,0);
        const Vector3d &pxn = xs.back().block<3,1>(0,3);
        // Initialize landmark position from estimated state
        p.tail<3>() = Rn*zi.second + pxn;
    } 
    zs.back().push_back(zi);
  }
}

} // namespace gcop
#endif
