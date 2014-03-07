#ifndef GCOP_BA_H
#define GCOP_BA_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include "pdmoc.h"
#include "bacost.h"
#include "kinbody2d.h"
#include "posegraph2d.h"

namespace gcop {
  
  using namespace Eigen;
  
  /**
   * Bundle adjustment based on optimal-control. This is a basic version
   * using a kinematic model in SE(2). It is possible to plug-in any dynamical model
   * without much overhead.
   */
  class Ba {
  public:
    
    /**
     * Bundle adjustment of pose-graph pg
     * @param pg pose-graph
     */
    Ba(Posegraph2d &pg);
    
    Posegraph2d &pg;   ///< pose graph
    
    Kinbody2d sys;     ///< system

    BaCost cost;       ///< cost function
    
    PDmoc<Matrix3d, 3, 3> pdmoc;  ///< parametric discrete-mechanics optimal control
  };
}

#endif
