#ifndef GCOP_BODY2DVIEW_H
#define GCOP_BODY2DVIEW_H

#include "body2d.h"
#include "systemview.h"

namespace gcop {
  using namespace Eigen;

  class Body2dView : public SystemView<pair<Matrix3d, Vector3d>, Vector3d> {
  public:

    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     */
    Body2dView(const Body2d &sys);
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    Body2dView(const Body2d &sys,
               vector<pair<Matrix3d, Vector3d> > *xs = 0,
               vector<Vector3d> *us = 0);
    
    virtual ~Body2dView();
    

    void Render(const pair<Matrix3d, Vector3d> *x,
                const Vector3d *u = 0);
   
    void Render(const vector<pair<Matrix3d, Vector3d> > *xs, 
                const vector<Vector3d> *us = 0,
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Body2d &sys;

    GLUquadricObj *qobj;

    double forceScale; ///< force arrow scaling coefficient
  };
}

#endif
