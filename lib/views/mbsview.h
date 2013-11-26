#ifndef GCOP_MBSVIEW_H
#define GCOP_MBSVIEW_H

#include "mbs.h"
#include "systemview.h"
#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#include "viewer.h"
#include "geom3dview.h"


namespace gcop {

#include "utils.h"

  using namespace Eigen;
  
    class MbsView : public SystemView<MbsState> {
  public:
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    MbsView(const Mbs &sys,
            vector<MbsState> *xs = 0);
    
    virtual ~MbsView();
    
    
    virtual void Render(const MbsState &x);
    
    void Render(const vector<MbsState> &xs, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Mbs &sys;
    GLUquadricObj *qobj;
    double dirSize;
    
    vector<Geom3dView*> geomViews;  
  };
}



#endif
