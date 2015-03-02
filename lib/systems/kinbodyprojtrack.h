#ifndef GCOP_KINBODYPROJTRACK_H
#define GCOP_KINBODYPROJTRACK_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include "kinbody3d.h"
#include "kinbody3dtrack.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  
  /**
   * Pose-track in 2D. This is the simplese possible definition assuming point (unprojected) features
   * e.g. from stereo/laser scans with uniform spherical covariance. The framework can
   * be easily extended to any general 2d/3d non-uniform or projected features as well.
   */
  class KinbodyProjTrack : public Kinbody3dTrack {
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
    KinbodyProjTrack(Kinbody3d &sys, int nf, double vd0, double t0, double tf,
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
  };
}

#endif
