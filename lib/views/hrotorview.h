#ifndef GCOP_HROTORVIEW_H
#define GCOP_HROTORVIEW_H

#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "body3dview.h"
#include "hrotorgeom3dview.h"
#include <Eigen/Dense>

namespace gcop {
  
  using namespace Eigen;
  
  class HrotorView : public Body3dView<4> {
  public:
  
  /**
   *  Create a particle view of trajectory traj
   * @param sys particle
   * @param xs trajectory
   */
  HrotorView(const Hrotor &sys,
             vector<pair<Matrix3d, Vector9d> > *xs = 0,
             vector<Vector4d> *us = 0);
  
  //  virtual ~HrotorView();  
  
  void Render(const pair<Matrix3d, Vector9d> *x,
              const Vector4d *u = 0);  

  HrotorGeom3dView view;
  };
}

#endif
