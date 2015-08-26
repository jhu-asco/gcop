#include "insview.h"
#include <iostream>
#include "so3.h"
#include "utils.h"

using namespace gcop;
using namespace Eigen;

InsView::InsView(const Ins &sys) : 
  SystemView<InsState, Vector6d>("Ins"), sys(sys)
{
  this->rgba[0] = .5;
  this->rgba[1] = .5;
  this->rgba[2] = .5;
  this->rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
  dirSize = -1;
}

InsView::InsView(const Ins &sys,
                 vector<InsState> *xs,
                 vector<Vector6d> *us) : 
  SystemView<InsState, Vector6d>("Ins", xs, us), sys(sys)
{
  this->rgba[0] = 0.5;
  this->rgba[1] = 0.5;
  this->rgba[2] = 0.5;
  this->rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
}


InsView::~InsView()
{
  free(qobj);
}

void InsView::Render(const InsState *x,
                     const Vector6d *u)
{
  //   glColor4f(1,0.5,0.5,0.5);  
  glPushMatrix();
  Transform(x->R, x->p);
  glScaled(.4, .2, .1); 
  glutSolidCube(1);
  glPopMatrix();
}


void InsView::Render(const vector<InsState> *xs, 
                     const vector<Vector6d> *us, 
                     bool rs, 
                     int is, int ie,
                     int dis, int dit,
                     bool dl)
{

  if (!xs->size()) {
    return;
  }
  
  // Viewer::SetColor(this->rgba[0], this->rgba[1], this->rgba[2], this->rgba[3]);
  glColor3f(this->rgba[0], this->rgba[1], this->rgba[2]);
  
  // set defaults
  if (is == -1)
    is = 0;
  if (ie == -1)
    ie = xs->size()-1;
  
  assert(is >= 0 && is <= xs->size()-1 && ie >= 0 && ie <= xs->size()-1);
  assert(is <= ie);
  
  glDisable(GL_LIGHTING);
  glLineWidth(this->lineWidth);
  glBegin(GL_LINE_STRIP);
  for (int i = is; i <= ie; i+=dit) {
    const InsState &x = (*xs)[i];
    glVertex3d(x.p[0], x.p[1], x.p[2]);
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

void InsView::Transform(const Matrix3d &R, const Vector3d &p)
{
  const SO3 &so3 = SO3::Instance();
  glTranslated(p[0], p[1], p[2]);
  
  Vector3d e;
  so3.log(e, R);
  double n = e.norm();
  if (n > so3.tol) {
    e = e/n;
    glRotated(RAD2DEG(n), e[0], e[1], e[2]);
  }
}
