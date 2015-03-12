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
    KinRccarPath(KinRccar &sys, bool odometry = false, bool forces = false);

    virtual void AddControl(const Vector2d &u, const double& v, double h);
    virtual void AddObservation(const vector< pair<int, Vector3d> > &z);

    KinRccar &sys;     ///< system

    bool odometry;
    bool extforce;         ///< is there a constant external force that need to be estimated
    bool forces;           ///< is there a constant external force that need to be estimated
    
    vector<double> ts;     ///< sequence of times (N+1 vector)

    vector<Matrix4d> xs;      ///< sequence of states x=(g,v) (N+1 vector)
    vector<Vector2d> us;   ///< sequence of inputs (N vector)
    vector<double> vs;   ///< sequence of odometry (N vector)
    
    vector<Matrix4d> xos;     ///< unoptimized trajectory
    vector<Vector2d> uos;     ///< unoptimized trajectory

    VectorXd p;            ///< observed landmark positions (2*m vector) and external force parameters (2 vector)    

    vector< vector<pair<int, Vector3d> > > zs;  ///< vector of vectors of landmark observations

    double cp;             ///< noise covariance of poses (assume spherical)
    double cv;           ///< noise covariance of odometry/gyro (assume diagonal)
    Vector2d cw;           ///< noise covariance of control forces (assume diagonal)

  };

KinRccarPath::KinRccarPath(KinRccar &sys, bool odometry, bool forces) : 
  sys(sys), odometry(odometry), forces(forces), extforce(false),
  p(extforce*3)
{
  Matrix4d x = MatrixXd::Identity(4,4);
  xs.push_back(x);
  xos.push_back(x);
  ts.push_back(0);

  cw = Vector2d(.05, pow(1.5*M_PI/180.,2));
  cv = .05;
  cp = .01;
}

void KinRccarPath::AddControl(const Vector2d &u, const double& v, double h)
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

void KinRccarPath::AddObservation(const vector< pair<int, Vector3d> > &z)
{
  zs.back() = z;
  // TODO: initialize p based on z from new lm
}

} // namespace gcop
#endif
