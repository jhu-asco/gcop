#ifndef GCOP_POSEGRAPH2D_H
#define GCOP_POSEGRAPH2D_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  /**
   * Pose-graph in 2D. This is the simplese possible definition assuming point (unprojected) features
   * e.g. from stereo/laser scans with uniform spherical covariance. The framework can
   * be easily extended to any general 2d/3d non-uniform or projected features as well.
   */
  class Posegraph2d {
  public:
    
    /**
     * Pose-landmark graph 
     * @param N number of trajectory segments
     * @param nf number of features (landmarks)
     */
    Posegraph2d(int N, int nf);    

    /**
     * Given a sequence of poses gs, compute the optimal feature locations p
     * @param p a vector of feature locations
     * @param gs a given vector of poses in SE(2)
     */
    void Optp(VectorXd &p, const vector<Matrix3d> &gs);

    /**
     * Create true and noisy synthetic pose graphs
     * @param pgt true pose graph
     * @param pgn noisy pose graph
     * @param tf total trajectory time
     */
    static void Synthesize(Posegraph2d &pgt, Posegraph2d &pgn, double tf);

    /**
     * Create true and noisy synthetic pose graphs
     * @param pgt true pose graph
     * @param pgn noisy pose graph
     * @param tf total trajectory time
     */
    static void Synthesize2(Posegraph2d &pgt, Posegraph2d &pgn, double tf);

    
    vector<double> ts;     ///< sequence of times (N+1 vector)

    vector<Matrix3d> gs;   ///< sequence of poses (N+1 vector)

    vector<Vector3d> us;   ///< sequence of inputs (N vector)

    VectorXd p;            ///< feature positions (2*nf vector)

    vector< vector<pair<int, Vector2d> > > Is;  ///< N+1-vector of vectors of visible pose-to-feature indices

    vector< vector<pair<int, Vector2d> > > Js;  ///< nf-vector of vectors of feature-to-pose indices indices
        
    double cp;             ///< noise covariance of poses (assume spherical)
  };
}

#endif
