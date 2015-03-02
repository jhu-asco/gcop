#ifndef GCOP_KINBODY3DTRACK_H
#define GCOP_KINBODY3DTRACK_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include "kinbody3d.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  
  /**
   * Pose-track in 2D. This is the simplese possible definition assuming point (unprojected) features
   * e.g. from stereo/laser scans with uniform spherical covariance. The framework can
   * be easily extended to any general 2d/3d non-uniform or projected features as well.
   */
  class Kinbody3dTrack {
  public:
    
    /**
     * Pose graph (using simulated features around a circular track)
     * @param sys body3d system
     * @param nf number of features (landmarks)
     * @param t0 start time
     * @param tf end time
     * @param r radius
     * @param odometry is odometry available
     * @param extforce treat constant external force in x-y as a parameter
     * @param forces use uncertain forces in the cost function
     */
    Kinbody3dTrack(Kinbody3d &sys, int nf, double vd0, double t0, double tf,
                double r = 25,
                bool extforce = false,
                bool forces = false);

    /**
     * Given a sequence of poses gs, compute the optimal feature locations p
     * @param p a vector of feature locations
     * @param gs a given vector of poses in SE(2)
     */
    virtual void Optp(VectorXd &p, const vector<Matrix4d> &xs);

    virtual void Get(Matrix4d &x, double vd, double t) const;

    virtual void MakeTrue();

    virtual void Add(const Vector6d &u, const Matrix4d &x, double h);

    /**
     * Add a new state/control to estimation vector
     * @param u commanded control 
     * @param x true state (to generate measurements)
     * @param h time-step 
     */
    virtual void Add2(const Vector6d &u, const Matrix4d &x, double h);

    Kinbody3d &sys;     ///< system

    double t0;       ///< initial time around track
    double tf;       ///< final time around track
    double r;       ///< track radius
    double w;      ///< track width
    double h;      ///< track height
    double dmax;    ///< sensor radius

    bool extforce;         ///< is there a constant external force that need to be estimated

    bool forces;           ///< is there a constant external force that need to be estimated
    
    vector<double> ts;     ///< sequence of times (N+1 vector)

    vector<Matrix4d> xs;      ///< sequence of states x=(g,v) (N+1 vector)

    vector<Vector6d> us;   ///< sequence of inputs (N vector)
    
    vector<Vector3d> ls;   ///< landmark positions
    vector<bool> observed; ///< was it observed

    vector<Matrix4d> xos;     ///< unoptimized trajectory
    vector<Vector6d> uos;     ///< unoptimized trajectory

    bool init; 
    VectorXd p;            ///< observed landmark positions (2*m vector) and external force parameters (2 vector)    
    vector<int> pis;       ///< observed landmark indices (m vector)

    double pr;     ///< feature radius

    vector< vector<pair<int, Vector3d> > > Is;  ///< N+1-vector of vectors of visible pose-to-feature indices

    vector< vector<pair<int, Vector3d> > > Js;  ///< nf-vector of vectors of feature-to-pose indices indices
    vector<int> cis;       ///< if observed then this should map into the corresponding value in p

    double cp;             ///< noise covariance of poses (assume spherical)

    Vector6d cv;           ///< noise covariance of odometry/gyro (assume diagonal)

    Vector6d cw;           ///< noise covariance of control forces (assume diagonal)

  };
}

#endif
