#ifndef GCOP_MBSCONTROLLER_H
#define GCOP_MBSCONTROLLER_H

#include "controller.h"
#include "mbs.h"

namespace gcop {

  using namespace std;
  using namespace Eigen;

  class MbsController : public Controller< MbsState, VectorXd> {
  public:

    /**
     * MBS PD controller
     *
     * @param sys multi-body system
     * @param xd desired state (optional, set to origin by default)
     * @param ad desired acceleration (optional, set to zero by default)
     */
    MbsController(const Mbs &sys,
                  MbsState *xd = 0, 
                  VectorXd *ad = 0);

    virtual void Set(VectorXd &u, double t, const MbsState &x);

    virtual ~MbsController();

    const Mbs &sys; ///< multi-body system

    MbsState *xd;   ///< desired state (origin by default)
    VectorXd *ad;   ///< desired acceleration (zero by default)

    VectorXd Kp;    ///< proportional terms  (ones by default)
    VectorXd Kd;    ///< derivative  terms  (ones by default)
    
  };
};

#endif
