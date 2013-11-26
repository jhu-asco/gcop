#ifndef GCOP_AIRMVIEW_H
#define GCOP_AIRMVIEW_H

#include "airm.h"
#include "cylinderview.h"
#include "mbsview.h"
#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#include "viewer.h"
#include "so3.h"
#include "hrotorgeom3dview.h"
#include "boxview.h"

namespace gcop {
  
  using namespace Eigen;
  
  class AirmView : public MbsView {
  public:
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    AirmView(const Airm &sys,
             vector<MbsState> *xs);
    
    HrotorGeom3dView hgv;
    //    BoxView views[2];
    CylinderView views[2];
  };
}

#endif
