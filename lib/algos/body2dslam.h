#ifndef GCOP_BODY2DSLAM_H
#define GCOP_BODY2DSLAM_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include "pddp.h"
#include "body2dslamcost.h"
#include "body2dgraph.h"

namespace gcop {
  
  using namespace Eigen;
  
  /**
   * Bundle adjustment based on optimal-control. This is a basic version
   * using a kinematic model in SE(2). It is possible to plug-in any dynamical model
   * without much overhead.
   */
  class Body2dSlam {
  public:
    
    /**
     * SLAM, global map iterative optimization
     * @param pg pose-graph
     */
    Body2dSlam(Body2dGraph &pg);
    
    Body2dGraph &pg;   ///< pose graph
    
    Body2dSlamCost cost;       ///< cost function
    
    PDdp<Body2dState, 6, 3> *pddp;  ///< parametric discrete-mechanics optimal control
  };
}

#endif
