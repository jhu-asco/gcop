#ifndef GCOP_GCARVIEW_H
#define GCOP_GCARVIEW_H

#include "gcar.h"
#include "systemview.h"

namespace gcop {
  using namespace Eigen;

  class GcarView : public SystemView<GcarState, Vector2d> {
  public:
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    GcarView(const Gcar &sys,
                  vector<GcarState> *xs = 0,
                  vector<Vector2d> *us = 0);

    virtual ~GcarView();
    
    
    void Render(const GcarState *x,
                const Vector2d *u = 0);
    
    void Render(const vector<GcarState> *xs, 
                const vector<Vector2d> *us, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Gcar &sys;
    GLUquadricObj *qobj;
  };
}

#endif
