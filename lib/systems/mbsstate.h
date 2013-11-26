#ifndef GCOP_MBSSTATE_H
#define GCOP_MBSSTATE_H

#include <Eigen/Dense>
#include "se3.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  // state dimension for nb-body system

  class MbsState {
  public:
  MbsState(int nb = 1) : 
    gs(nb), vs(nb), dgs(nb-1),
      r(nb-1), dr(nb-1) {
    }

  MbsState(const MbsState &x) : 
    gs(x.gs), vs(x.vs), dgs(x.dgs),
      r(x.r), dr(x.dr) {
    }

    
    vector<Matrix4d> gs;        ///< configurations
    vector<Vector6d> vs;        ///< body-vixed velocities
    vector<Matrix4d> dgs;       ///< relative xforms from b/n bodies
    
    VectorXd r;  ///< joint angles
    VectorXd dr; ///< joint velocities
  };
}

#endif
