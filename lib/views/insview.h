#ifndef GCOP_INSVIEW_H
#define GCOP_INSVIEW_H

#include "ins.h"
#include "systemview.h"
#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

namespace gcop {

  using namespace Eigen;
  
  class InsView : public SystemView<InsState, Vector6d> {
    
  public:
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     */    
    InsView(const Ins& sys);
    
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    InsView(const Ins &sys,
            vector<InsState> *xs,
            vector<Vector6d> *us = 0);
    
    virtual ~InsView();
    
  
    virtual void Render(const InsState *x,
                        const Vector6d *u = 0);
    
    void Render(const vector<InsState> *xs, 
                const vector<Vector6d> *us = 0, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Ins &sys;
    GLUquadricObj *qobj;
    double dirSize;
    
    static void Transform(const Matrix3d &R, const Vector3d &p);
  };
}

#endif
