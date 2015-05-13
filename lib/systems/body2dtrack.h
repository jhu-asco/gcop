#ifndef GCOP_BODY2DTRACK_H
#define GCOP_BODY2DTRACK_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include "body2d.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  
  /**
   * Pose-track in 2D. This is the simplese possible definition assuming point (unprojected) features
   * e.g. from stereo/laser scans with uniform spherical covariance. The framework can
   * be easily extended to any general 2d/3d non-uniform or projected features as well.
   */
  class Body2dTrack {
  public:
    
    /**
     * Pose graph (using simulated features around a circular track)
     * @param sys body2d system
     * @param nf number of features (landmarks)
     * @param t0 start time
     * @param tf end time
     * @param r radius
     * @param odometry is odometry available
     * @param extforce treat constant external force in x-y as a parameter
     * @param forces use uncertain forces in the cost function
     */
    Body2dTrack(Body2d<> &sys, int nf, double t0, double tf,
                double r = 25,
                bool odometry = true,
                bool extforce = false,
                bool forces = false);

    /**
     * Given a sequence of poses gs, compute the optimal feature locations p
     * @param p a vector of feature locations
     * @param gs a given vector of poses in SE(2)
     */
    void Optp(VectorXd &p, const vector<Body2dState> &xs);

    void Get(Body2dState &x, double vd, double t) const;

    void MakeTrue();

    void Add(const Vector3d &u, const Body2dState &x, double h);

    /**
     * Add a new state/control to estimation vector
     * @param u commanded control 
     * @param x true state (to generate measurements)
     * @param h time-step 
     */
    void Add2(const Vector3d &u, const Body2dState &x, double h);

    Body2d<> &sys;     ///< system

    double t0;       ///< initial time around track
    double tf;       ///< final time around track
    double r;       ///< track radius
    double w;      ///< track width
    double dmax;    ///< sensor radius

    bool odometry;         ///< is there odometry available (this means that the velocity part of the state will be regarded as available measurements)

    bool extforce;         ///< is there a constant external force that need to be estimated

    bool forces;           ///< is there a constant external force that need to be estimated
    
    vector<double> ts;     ///< sequence of times (N+1 vector)

    vector<Body2dState> xs;      ///< sequence of states x=(g,v) (N+1 vector)

    vector<Vector3d> us;   ///< sequence of inputs (N vector)
    
    vector<Vector2d> ls;   ///< landmark positions
    vector<bool> observed; ///< was it observed

    vector<Body2dState> xos;     ///< unoptimized trajectory
    vector<Vector3d> uos;     ///< unoptimized trajectory

    bool init; 
    VectorXd p;            ///< observed landmark positions (2*m vector) and external force parameters (2 vector)    
    vector<int> pis;       ///< observed landmark indices (m vector)

    double pr;     ///< feature radius

    vector< vector<pair<int, Vector2d> > > Is;  ///< N+1-vector of vectors of visible pose-to-feature indices

    vector< vector<pair<int, Vector2d> > > Js;  ///< nf-vector of vectors of feature-to-pose indices indices
    vector<int> cis;       ///< if observed then this should map into the corresponding value in p

        
    vector<Vector3d> vs;   ///< odometry/gyro velocity measurements

    double cp;             ///< noise covariance of poses (assume spherical)

    Vector3d cv;           ///< noise covariance of odometry/gyro (assume diagonal)

    Vector3d cw;           ///< noise covariance of control forces (assume diagonal)

  };
}

#endif
