#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "gcarview.h"
#include "viewer.h"
#include "se2.h"
#include <math.h>

using namespace gcop;
using namespace Eigen;

GcarView::GcarView(const Gcar &sys, 
                   vector< pair<Matrix3d, double> > *xs,
                   vector<Vector2d> *us) : 
  SystemView("Gcar", xs, us), sys(sys)
{
  rgba[0] = 0.5;
  rgba[1] = 0.5;
  rgba[2] = 0.5;
  rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
}


GcarView::~GcarView()
{
  free(qobj);
}


void GcarView::Render(const pair<Matrix3d, double> *x,
                           const Vector2d *u)
{
  Vector3d q;
  SE2::Instance().g2q(q, x->first);
  const double &theta = q[0];
  const double &px = q[1];
  const double &py = q[2];
  //  const double &w = x->second[0];
  const double &v = x->second;

  double phi = u ? atan((*u)[1]) : 0;
  double d = sys.l;

  glPushMatrix();
  glTranslated(px, py, d/8);
  glRotated(RAD2DEG(theta), 0,0,1);

  glTranslated(d/2, 0, 0);

  glPushMatrix();
  glScaled(1.2*d, .6*d, .2*d);
  glBegin(GL_LINE_LOOP);
  glVertex2d(-.5,-.5);
  glVertex2d(.5,-.5);
  glVertex2d(.5,.5);
  glVertex2d(-.5,.5);
  glEnd();
  //  glutSolidCube(1);
  glPopMatrix();
  
  glPushMatrix();
  glTranslated(d/2,d/4-d/16,0);
  glRotated(RAD2DEG(phi),0,0,1);
  glRotated(90,1,0,0);
  glTranslated(0,0,-d/16);
  gluCylinder(qobj, d/8, d/8, d/8, 10, 10);
  glPopMatrix();

  glPushMatrix();
  glTranslated(d/2,-d/4+d/16,0);
  glRotated(RAD2DEG(phi),0,0,1);
  glRotated(90,1,0,0);
  glTranslated(0,0,-d/16);
  gluCylinder(qobj, d/8, d/8, d/8, 10, 10);
  glPopMatrix();

  glPushMatrix();
  glTranslated(-d/2,d/4-d/16,0);
  glRotated(90,1,0,0);
  glTranslated(0,0,-d/16);
  gluCylinder(qobj, d/8, d/8, d/8, 10, 10);
  glPopMatrix();

  glPushMatrix();
  glTranslated(-d/2,-d/4+d/16,0);
  glRotated(90,1,0,0);
  glTranslated(0,0,-d/16);
  gluCylinder(qobj, d/8, d/8, d/8, 10, 10);
  glPopMatrix();

  glPopMatrix();
}


void GcarView::Render(const vector<pair<Matrix3d, double> > *xs,
                           const vector<Vector2d> *us,
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
    const pair<Matrix3d, double> &x = (*xs)[i];
    glVertex3d(x.first(0,2), x.first(1,2), 0);
  }
  glEnd();
  glLineWidth(1);
  glEnable(GL_LIGHTING);
  
  if (xs->size())
    if (us && us->size())
      Render(&(*xs)[0], &(*us)[0]);
    else
      Render(&(*xs)[0]);

  if (rs) {    
    for (int i = 1; i < xs->size(); i+=dis) {
      if (i < us->size())
        Render(&(*xs)[i], &(*us)[i]);
      else
        Render(&(*xs)[i]);
    }    
  }
  
  if (dl)
    Render(&xs->back());
}
