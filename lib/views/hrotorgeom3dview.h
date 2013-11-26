#ifndef GCOP_HROTORGEOM3DVIEW_H
#define GCOP_HROTORGEOM3DVIEW_H

#include "geom3dview.h"
#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "hrotor.h"


namespace gcop {
  
  using namespace Eigen;
  
  class HrotorGeom3dView : public Geom3dView {
  public:    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     */    
    
    HrotorGeom3dView(const Hrotor &sys,
                     Matrix4d *g = 0);
    
    virtual ~HrotorGeom3dView();
    
    void RenderGeom();
    
    const Hrotor &sys;

    GLUquadricObj *qobj;

  };
}

#endif
