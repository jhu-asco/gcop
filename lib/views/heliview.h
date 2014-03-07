#ifndef GCOP_HELIVIEW_H
#define GCOP_HELIVIEW_H

#include "heli.h"
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
  
  class HeliView : public Body3dView<4> {
  public:
  
  /**
   *  Create a particle view of trajectory traj
   * @param sys particle
   */
  
  HeliView(const Heli &sys);
  

  virtual ~HeliView();
  
  /**
   *  Create a particle view of trajectory traj
   * @param sys particle
   * @param xs trajectory
   */
  HeliView(const Heli &sys,
           vector<pair<Matrix3d, Vector9d> > *xs,
           vector<Vector4d> *us = 0);
  
  //  virtual ~HeliView();  
  
  void Render(const pair<Matrix3d, Vector9d> *x,
              const Vector4d *u = 0);
  
  const Heli &sys;

    GLUquadricObj *body;
    GLUquadricObj *tail;
    GLUquadricObj* tprop[4];
    GLUquadricObj *rprop;


  };
}
#endif
