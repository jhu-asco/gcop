#ifndef GCOP_KINBODY2DVIEW_H
#define GCOP_KINBODY2DVIEW_H

#include "kinbody2d.h"
#include "systemview.h"

namespace gcop {
  using namespace Eigen;

  class Kinbody2dView : public SystemView<Matrix3d> {
  public:

    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     */
    Kinbody2dView(const Kinbody2d &sys);
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    Kinbody2dView(const Kinbody2d &sys,
                  vector<Matrix3d> *xs);

    virtual ~Kinbody2dView();
       

    void Render(const Matrix3d &x);
    
    void Render(const vector<Matrix3d> &xs, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Kinbody2d &sys;

    GLUquadricObj *qobj;
  };
}

#endif
