#ifndef GCOP_UNICYCLEVIEW_H
#define GCOP_UNICYCLEVIEW_H

#include "unicycle.h"
#include "systemview.h"

namespace gcop {
  using namespace Eigen;

  class UnicycleView : public SystemView<Vector5d, Vector2d> {
  public:
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    UnicycleView(const Unicycle &sys,
                 vector<Vector5d> *xs = 0,
                 vector<Vector2d> *us = 0);

    virtual ~UnicycleView();
       

    void Render(const Vector5d *x = 0,
                const Vector2d *u = 0);
    
    void Render(const vector<Vector5d> *xs, 
                const vector<Vector2d> *us = 0, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Unicycle &sys;
    GLUquadricObj *qobj;
  };
}

#endif
