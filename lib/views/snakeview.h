#ifndef GCOP_SNAKEVIEW_H
#define GCOP_SNAKEVIEW_H

#include "snake.h"
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
  
  class SnakeView : public MbsView {
  public:
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    SnakeView(const Snake &sys,
              vector<MbsState> *xs);
    
    
    BoxView views[3];
  };
}
#endif
