#ifndef GCOP_POINT3DVIEW_H
#define GCOP_POINT3DVIEW_H

#include "point3d.h"
#include "systemview.h"

namespace gcop {
  using namespace Eigen;

  class Point3dView : public SystemView<Point3dState, Vector3d> {
  public:
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    Point3dView(const Point3d &sys,
                vector<Point3dState> *xs = 0,
                vector<Vector3d> *us = 0);
    
    virtual ~Point3dView();
       
    
    void Render(const Point3dState *x,
                const Vector3d *u = 0);
    
    void Render(const vector<Point3dState> *xs,
                const vector<Vector3d> *us = 0,
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Point3d &sys;

    double r; ///< radius
  };
}

#endif
