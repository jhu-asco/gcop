#ifndef GCOP_AIRBOTVIEW_H
#define GCOP_AIRBOTVIEW_H

#include "airbot.h"
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
  
  class AirbotView : public MbsView {
  public:
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    AirbotView(const Airbot &sys,
                vector<MbsState> *xs);
    
    HrotorGeom3dView hgv;
    //    BoxView views[2];
    CylinderView views[6];
  };
}

#endif
