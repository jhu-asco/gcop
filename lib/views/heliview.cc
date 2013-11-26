#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "heliview.h"
#include "viewer.h"
#include "so3.h"
#include "utils.h"

using namespace gcop;
using namespace Eigen;

HeliView::HeliView(const Heli &sys) : 
  Body3dView(sys), sys(sys)
{
}

HeliView::HeliView(const Heli &sys, 
                   vector< pair<Matrix3d, Vector9d> > *xs) : 
  Body3dView(sys, xs), sys(sys)
{
  body = gluNewQuadric();
  tail = gluNewQuadric();
  for (int i=0;i<4;++i)
    tprop[i] = gluNewQuadric();
  rprop = gluNewQuadric();
}

HeliView::~HeliView()
{
  free(rprop);
  for (int i = 0; i < 4; ++i)
    free(tprop[i]);
  free(tail);
  free(body);
}

void HeliView::Render(const pair<Matrix3d, Vector9d> &x)
{
  //   glColor4f(1,0.5,0.5,0.5);

  glPushMatrix();
  Transform(x.first, x.second.head<3>());
  
  Viewer::SetColor(0,0.5,0.5,0);

  //  glColor3f(0,.5,.5);

  glutSolidSphere(.9*sys.rt, 10, 10);

  //  glColor3f(0,.5,.5);
  Viewer::SetColor(0, 0.5, 0.5, 0);

  glPushMatrix();
  glRotated(-90, 0, 1, 0);
  gluCylinder(tail, sys.rt/4, sys.rt/8, sys.rb, 10, 10);
  glPopMatrix();  
  
  //  glColor3f(0,0,.5);
  Viewer::SetColor(0, 0, 0.5, 0);

  glPushMatrix();
  glTranslated(-sys.rb, 0, 0);
  //  glBegin(GL_LINE);
  //glVertex3d(0,0,0);
  //glVertex3d(0,-heli.us[4*i+1]/10,0);
  //glEnd();
  glRotated(-90, 1, 0, 0);
  //gluCylinder(rprop, sys.rt/2, sys.rt/2, .01, 10, 10);
  gluDisk(rprop, 0, sys.rt/2, 10, 10);
  glPopMatrix();
  

  //  glColor4f(0,0,.5, .2);
  Viewer::SetColor(0, 0, 0.5, 0.2);

  glPushMatrix();
  glTranslated(0,0,sys.rt);

  //  glRotated(RAD2DEG(ar), 1, 0, 0);
  //  glRotated(RAD2DEG(ap), 0, 1, 0);

  int i0 = (int)round(rand()/(double)RAND_MAX);
  for (int j=0;j<4;++j)
    gluPartialDisk(rprop, 0, sys.rb, 10, 10, i0*45 + 22.5+j*90, 45);
  //gluCylinder(tprop, sys.rb, sys.rb, .05, 10, 10);
  glPopMatrix();

  glPopMatrix();
}
