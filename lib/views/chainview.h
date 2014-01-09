#ifndef GCOP_CHAINVIEW_H
#define GCOP_CHAINVIEW_H

#include "chain.h"
#include "boxview.h"
#include "mbsview.h"
#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#include "viewer.h"
#include "so3.h"

namespace gcop {
  
  using namespace Eigen;
  
  class ChainView : public MbsView {
  public:
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    ChainView(const Chain &sys,
              vector<MbsState> *xs);
        
    virtual ~ChainView();

    BoxView **views;
  };
}
#endif
