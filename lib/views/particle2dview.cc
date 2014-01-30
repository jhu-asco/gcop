#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "particle2dview.h"
#include "utils.h"
#include "viewer.h"

using namespace gcop;
using namespace Eigen;

Particle2dView::Particle2dView(const Particle2d &sys, 
                               vector<Vector4d> *xs,
                               vector<Vector2d> *us) : 
  SystemView("Particle", xs, us), sys(sys)
{
  rgba[0] = 0.5;
  rgba[1] = 0.5;
  rgba[2] = 0.5;
  rgba[3] = 0;
  this->lineWidth = 2;
}


Particle2dView::~Particle2dView()
{
}


void Particle2dView::Render(const Vector4d *x, const Vector2d *u)
{
  glPushMatrix();
  glTranslated((*x)[0], (*x)[1],  0);  
  glutSolidSphere(sys.r, 10, 10);
  glPopMatrix();  
}

void Particle2dView::Render(const vector<Vector4d> *xs, 
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
    const Vector4d &x = (*xs)[i];
    glVertex3d(x[0], x[1], 0);
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
