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
                       vector< pair<Matrix3d, Vector9d> > *xs) : 
  Body3dView(sys, xs), view(sys)
{
}

void HrotorView::Render(const pair<Matrix3d, Vector9d> &x)
{
  //   glColor4f(1,0.5,0.5,0.5);
  
  glPushMatrix();
  Transform(x.first, x.second.head<3>());
  view.RenderGeom();
  glPopMatrix();
}
