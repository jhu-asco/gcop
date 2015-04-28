#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "point3dview.h"
#include "utils.h"
#include "viewer.h"

using namespace gcop;
using namespace Eigen;

Point3dView::Point3dView(const Point3d &sys, 
                         vector<Point3dState> *xs,
                         vector<Vector3d> *us) : 
  SystemView("Particle", xs, us), sys(sys)
{
  rgba[0] = 0.5;
  rgba[1] = 0.5;
  rgba[2] = 0.5;
  rgba[3] = 0;
  this->lineWidth = 2;
  this->r = .1;
}


Point3dView::~Point3dView()
{
}


void Point3dView::Render(const Point3dState *x, const Vector3d *u)
{
  assert(x);
  glPushMatrix();
  glTranslated(x->q[0], x->q[1], x->q[2]);
  glutSolidSphere(this->r, 10, 10);
  glPopMatrix();  
}

void Point3dView::Render(const vector<Point3dState> *xs, 
                            const vector<Vector3d> *us, 
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
    ie = xs->size()-1;

  assert(is >= 0 && is <= xs->size()-1 && ie >= 0 && ie <= xs->size()-1);
  assert(is <= ie);

  glDisable(GL_LIGHTING);
  glLineWidth(lineWidth);
  glBegin(GL_LINE_STRIP);
  for (int i = is; i <= ie; i+=dit) {
    const Point3dState &x = (*xs)[i];
    glVertex3d(x.q[0], x.q[1], x.q[2]);
  }
  glEnd();
  glLineWidth(1);
  glEnable(GL_LIGHTING);
  
  if (rs) {
    for (int i = 0; i < xs->size(); i+=dis) {
      Render(&(*xs)[i]);
    }
  }

  if (dl)
    Render(&xs->back());
}
