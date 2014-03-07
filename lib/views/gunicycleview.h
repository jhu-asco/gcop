#ifndef GCOP_GUNICYCLEVIEW_H
#define GCOP_GUNICYCLEVIEW_H

#include "gunicycle.h"
#include "systemview.h"

namespace gcop {
  using namespace Eigen;

  class GunicycleView : public SystemView<pair<Matrix3d, Vector2d>, Vector2d> {
  public:
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    GunicycleView(const Gunicycle &sys,
                  vector<pair<Matrix3d, Vector2d> > *xs = 0,
                  vector<Vector2d> *us = 0);

    virtual ~GunicycleView();
    
    
    void Render(const pair<Matrix3d, Vector2d> *x,
                const Vector2d *u = 0);
    
    void Render(const vector<pair<Matrix3d, Vector2d> > *xs, 
                const vector<Vector2d> *us, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Gunicycle &sys;
    GLUquadricObj *qobj;
  };
}

#endif
