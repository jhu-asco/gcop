#ifndef GCOP_GCARVIEW_H
#define GCOP_GCARVIEW_H

#include "gcar.h"
#include "systemview.h"

namespace gcop {
  using namespace Eigen;

  class GcarView : public SystemView<pair<Matrix3d, double>, Vector2d> {
  public:
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    GcarView(const Gcar &sys,
                  vector<pair<Matrix3d, double> > *xs = 0,
                  vector<Vector2d> *us = 0);

    virtual ~GcarView();
    
    
    void Render(const pair<Matrix3d, double> *x,
                const Vector2d *u = 0);
    
    void Render(const vector<pair<Matrix3d, double> > *xs, 
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
