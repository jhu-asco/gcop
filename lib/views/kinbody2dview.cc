#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "kinbody2dview.h"
#include "viewer.h"
#include "se2.h"

using namespace gcop;
using namespace Eigen;

Kinbody2dView::Kinbody2dView(const Kinbody2d &sys) : 
  SystemView("Kinbody2d"), sys(sys)
{
  rgba[0] = .5;
  rgba[1] = .5;
  rgba[2] = .5;
  rgba[3] = 0.5;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
}

Kinbody2dView::Kinbody2dView(const Kinbody2d &sys, 
                             vector<Matrix3d> *xs) : 
  SystemView("Kinbody2d", xs), sys(sys)
{
  rgba[0] = 0.5;
  rgba[1] = 0.5;
  rgba[2] = 0.5;
  rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
}


Kinbody2dView::~Kinbody2dView()
{
  free(qobj);
}


void Kinbody2dView::Render(const Matrix3d &x)
{
  Vector3d q;
  SE2::Instance().g2q(q, x);

  Viewer::SetColor(rgba[0], rgba[1], rgba[2], rgba[3]);
  glPushMatrix();
  glTranslated(q[1], q[2], 0);
  glRotated(RAD2DEG(q[0]), 0,0,1);
  glScaled(sys.d[0], sys.d[1], sys.d[1]/4); 
  glutSolidCube(1);
  glPopMatrix();
}


void Kinbody2dView::Render(const vector<Matrix3d> &xs, 
                           bool rs, 
                           int is, int ie,
                           int dis, int dit,
                           bool dl)
{
  Viewer::SetColor(rgba[0], rgba[1], rgba[2], rgba[3]);
  //  glColor4f(

  // set defaults
  if (is == -1)
    is = 0;
  if (ie == -1)
    ie = xs.size()-1;

  assert(is >= 0 && is <= xs.size()-1 && ie >= 0 && ie <= xs.size()-1);
  assert(is <= ie);

  glDisable(GL_LIGHTING);

  glLineWidth(lineWidth);
  glBegin(GL_LINE_STRIP);
  for (int i = is; i <= ie; i+=dit) {
    const Matrix3d &x = xs[i];
    glVertex3d(x(0,2), x(1,2), 0);
  }
  glEnd();
  glLineWidth(1);

    glEnable(GL_LIGHTING);
  
  if (rs) {
    for (int i = 0; i < xs.size(); i+=dis) {
      Render(xs[i]);
    }
  }

  if (dl)
    Render(xs.back());
}
