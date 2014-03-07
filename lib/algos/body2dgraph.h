#ifndef GCOP_BODY2DGRAPH_H
#define GCOP_BODY2DGRAPH_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include "body2d.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  /**
   * Pose-graph in 2D. This is the simplese possible definition assuming point (unprojected) features
   * e.g. from stereo/laser scans with uniform spherical covariance. The framework can
   * be easily extended to any general 2d/3d non-uniform or projected features as well.
   */
  class Body2dGraph {
  public:
    
    /**
     * Pose-landmark graph 
     * @param sys body2d system
     * @param N number of trajectory segments
     * @param nf number of features (landmarks)
     * @param odometry is odometry available
     * @param extforce treat constant external force in x-y as a parameter
     * @param forces use uncertain forces in the cost function
     */
    Body2dGraph(Body2d &sys, int N, int nf, 
                bool odometry = true, 
                bool extforce = false,
                bool forces = false);

    /**
     * Given a sequence of poses gs, compute the optimal feature locations p
     * @param p a vector of feature locations
     * @param gs a given vector of poses in SE(2)
     */
    void Optp(VectorXd &p, const vector<M3V3d> &xs);

    /**
     * Create true and noisy synthetic pose graphs
     * @param pgt true pose graph
     * @param pgn noisy pose graph
     * @param tf total trajectory time
     */
    static void Synthesize(Body2dGraph &pgt, Body2dGraph &pgn, double tf);

    /**
     * Create true and noisy synthetic pose graphs
     * @param pgt true pose graph
     * @param pgn noisy pose graph
     * @param tf total trajectory time
     */
    static void Synthesize2(Body2dGraph &pgt, Body2dGraph &pgn, double tf);

    /**
     * Create true and noisy synthetic pose graphs
     * @param pgt true pose graph
     * @param pgn noisy pose graph
     * @param tf total trajectory time
     */
    static void Synthesize3(Body2dGraph &pgt, Body2dGraph &pgn, double tf);

    Body2d &sys;     ///< system

    bool odometry;         ///< is there odometry available (this means that the velocity part of the state will be regarded as available measurements)
    
    bool extforce;         ///< is there a constant external force that need to be estimated

    bool forces;         ///< optimize over forces

    Vector3d  fext;      ///< unknown external force 
    
    vector<double> ts;     ///< sequence of times (N+1 vector)

    vector<M3V3d>  xs;      ///< sequence of states x=(g,v) (N+1 vector)

    vector<Vector3d> us;   ///< sequence of inputs (N vector)

    VectorXd p;            ///< feature positions (2*nf vector) and external force parameters (2 vector)

    vector< vector<pair<int, Vector2d> > > Is;  ///< N+1-vector of vectors of visible pose-to-feature indices

    vector< vector<pair<int, Vector2d> > > Js;  ///< nf-vector of vectors of feature-to-pose indices indices
        
    vector<Vector3d> vs;   ///< odometry/gyro velocity measurements

    double cp;             ///< noise covariance of poses (assume spherical)

    Vector3d cv;           ///< noise covariance of odometry/gyro (assume diagonal)

    Vector3d cw;           ///< noise covariance of control forces (assume diagonal)

  };
}

#endif
