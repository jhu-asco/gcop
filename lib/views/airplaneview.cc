#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "airplaneview.h"
#include "viewer.h"
#include "se2.h"

using namespace gcop;
using namespace Eigen;

AirplaneView::AirplaneView(const Gunicycle &sys, 
                           vector< pair<Matrix3d, Vector2d> > *xs) : 
  SystemView("Airplane", xs), sys(sys)
{
  rgba[0] = 0.5;
  rgba[1] = 0.5;
  rgba[2] = 0.5;
  rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
}


AirplaneView::~AirplaneView()
{
  free(qobj);
}


void AirplaneView::Render(const pair<Matrix3d, Vector2d> &x)
{
  Vector3d q;
  SE2::Instance().g2q(q, x.first);
  const double &theta = q[0];
  const double &px = q[1];
  const double &py = q[2];
  const double &w = x.second[0];
  const double &v = x.second[1];

  double phi = 0;
  double d = sys.dx;

  glPushMatrix();
  
  glTranslated(px, py, 0);
  glRotated(RAD2DEG(theta), 0,0,1);

  glScalef(.3, .3, .3);

  Viewer::SetColor(0, 0, 1, 0);

  glPushMatrix();
  glTranslated(1,0,0);
  glRotated(-90, 0, 1, 0);
  glScaled(.2,1,1);
  gluCylinder(qobj, .66, 1.5, 2, 10, 10);
  glPopMatrix();  
  
  glPopMatrix();

}


void AirplaneView::Render(const vector<pair<Matrix3d, Vector2d> > &xs, 
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
    const pair<Matrix3d, Vector2d> &x = xs[i];
    glVertex3d(x.first(0,2), x.first(1,2), 0);
  }
  glEnd();
  glLineWidth(1);
  glEnable(GL_LIGHTING);
  
  if (xs.size())
    Render(xs[0]);
  
  if (rs) {    
    for (int i = 1; i < xs.size(); i+=dis) {
      Render(xs[i]);
    }    
  }
  
  if (dl)
    Render(xs.back());
}
