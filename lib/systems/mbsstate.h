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
  MbsState(int nb = 1, bool fixed = false) : 
    gs(nb), vs(nb), dgs(nb-1),
      r(nb-1), dr(nb-1), ub(nb-1, false), lb(nb-1, false), 
      zu(VectorXd::Zero(nb-1)), zl(VectorXd::Zero(nb-1)), fixed(fixed) {
    }
    
  MbsState(const MbsState &x) : 
    gs(x.gs), vs(x.vs), dgs(x.dgs),
      r(x.r), dr(x.dr), ub(x.ub), lb(x.lb), 
      zl(x.zl), zu(x.zu), fixed(x.fixed) {
    }
    
    vector<Matrix4d> gs;        ///< configurations
    vector<Vector6d> vs;        ///< body-fixed velocities
    vector<Matrix4d> dgs;       ///< relative xforms from b/n bodies
    
    VectorXd r;  ///< joint angles
    VectorXd dr; ///< joint velocities
    vector<bool> ub; ///< at upper bound
    vector<bool> lb; ///< at lower bound

    VectorXd zl;  ///< lower bound spring
    VectorXd zu;  ///< upper bound spring

    bool fixed;  ///< fixed base?
  };
}

#endif
