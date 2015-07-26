#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "hrotorview.h"
#include "viewer.h"
#include "so3.h"
#include "utils.h"

using namespace gcop;
using namespace Eigen;

HrotorView::HrotorView(const Hrotor &sys, 
                       vector< Body3dState > *xs,
                       vector<Vector4d> *us) : 
  Body3dView(sys, xs, us), view(sys)
{
}

void HrotorView::Render(const Body3dState *x,
                        const Vector4d *u)
{
  //   glColor4f(1,0.5,0.5,0.5);
  
  glPushMatrix();
  Transform(x->R, x->p);
  view.RenderGeom();
  glPopMatrix();
}
