#ifndef GCOP_QROTORVIEW_H
#define GCOP_QROTORVIEW_H

#include "qrotor.h"
#include "body3dview.h"
#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#include "viewer.h"
#include "so3.h"


namespace gcop {
  
  using namespace Eigen;
  
  class QrotorView : public Body3dView<4> {
  public:

  /**
   *  Create a particle view of trajectory traj
   * @param sys particle
   * @param xs trajectory
   */
  QrotorView(const Qrotor &sys,
              vector<pair<Matrix3d, Vector9d> > *xs = 0);
  
  //  virtual ~QrotorView();  
  
  void Render(const pair<Matrix3d, Vector9d> &x);

  const Qrotor &sys;
  };
}
#endif
