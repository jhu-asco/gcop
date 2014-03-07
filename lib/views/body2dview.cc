#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "body2dview.h"
#include "viewer.h"
#include "se2.h"

using namespace gcop;
using namespace Eigen;

Body2dView::Body2dView(const Body2d &sys) : 
  SystemView("Body2d"), sys(sys)
{
  rgba[0] = .5;
  rgba[1] = .5;
  rgba[2] = .5;
  rgba[3] = 0.5;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
  forceScale = 1;
}

Body2dView::Body2dView(const Body2d &sys, 
                       vector< pair<Matrix3d, Vector3d> > *xs,
                       vector<Vector3d> *us) : 
  SystemView("Body2d", xs, us), sys(sys) 
{
  rgba[0] = 0.5;
  rgba[1] = 0.5;
  rgba[2] = 0.5;
  rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
  forceScale = 1;
}


Body2dView::~Body2dView()
{
  free(qobj);
}


static void RenderUuv(GLUquadricObj *qobj, double l, double c)
{

  glPushMatrix();

  glRotated(90, 0,1,0);
  glTranslated(0,0,-l);

  // main body
  //  glColor3f(.9,.9,.9);
  gluCylinder(qobj, 0.6*c, 0.6*c, 2*l, 10, 10);
  gluDisk(qobj, c/6, c/2, 10, 10);

  // tip
  //  glColor3f(0,.4,.4);
  glPushMatrix();
  glTranslated(0,0,2*l);
  glutSolidSphere(0.58*c, 10, 10);
  glPopMatrix();

  //  glColor3f(.4,.4,.4);
  /*
  glPushMatrix();
  glTranslated(0,0,c/3);
  glScaled(c*.33, c*1.6, c*.6);
  glutSolidSphere(1, 10, 10);
  glPopMatrix();
  */

  double r = c;

    glColor3f(0,0,0);
    
    glPushMatrix();
    
    glTranslated(0, 0, -.7*c);
    glScaled(.7,.7,.7);
    // motor
    glPushMatrix();
    glTranslated(0,0,r);
    glScaled(.5,.5,1);
    glutSolidSphere(r, 10, 10);      

    //    glColor4f(0, 0, .5, .8);
    glPopMatrix();

    // motor stand
    glPushMatrix();
    glTranslated(-r,0,r);
    glScaled(1,.1,.1);

    glColor3f(.1,.1,.1);
    glutSolidSphere(.5*c, 10, 10);

    glPopMatrix();

    gluCylinder(qobj, r, r, r, 10, 10);
    
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);      
    
    glPushMatrix();    
    glColor4f(0, 0, .5,.3);

    //    gluCylinder(qobj, fabs(s.u[i])*.03, fabs(s.u[i])*.2, fabs(s.u[i]), 10, 10);
    //      std::cout << "PUT" << std::endl;
      //    }
    glPopMatrix();
    glDisable(GL_BLEND);      

    //  glPushMatrix();
 
    //    glColor(.3,.3,.3);

    glTranslated(0,0,l/16);
    
    //  glRotated(RAD2DEG(ar), 1, 0, 0);
    //  glRotated(RAD2DEG(ap), 0, 1, 0);
    
    //    int i0 = (int)round(rand()/(double)RAND_MAX);
    int i0=0;
    for (int j = 0; j < 4; ++j)
      gluPartialDisk(qobj, 0, r, 10, 10, i0*45 + 22.5+j*90, 22.5);
    //gluCylinder(tprop, rb, rb, .05, 10, 10);
    
    glPopMatrix();    
    
  
  glPopMatrix();
}



void Body2dView::Render(const pair<Matrix3d, Vector3d> *x, const Vector3d *u)
{
  Vector3d q;
  SE2::Instance().g2q(q, x->first);

  glDisable(GL_LIGHTING);  
  //  Viewer::SetColor(rgba[0], rgba[1], rgba[2], rgba[3]);
  glPushMatrix();
  glTranslated(q[1], q[2], 0);
  glRotated(RAD2DEG(q[0]), 0,0,1);
  //  glScaled(sys.d[0], sys.d[1],  sys.d[1]/4);
  //  glutSolidSphere(1, 10, 10);
  //  glutSolidCube(1);


  glColor4d(rgba[0], rgba[1], rgba[2], rgba[3]);
  RenderUuv(qobj, 1,  .4);
  glColor4d(rgba[0], rgba[1], rgba[2], rgba[3]);
  if (renderForces && u) {
    Vector3d f(forceScale*(*u)[1], forceScale*(*u)[2], 0);
    Viewer::DrawArrow(f.data(), qobj);
  }

  glPopMatrix();
  glEnable(GL_LIGHTING);  

}


void Body2dView::Render(const vector<pair<Matrix3d, Vector3d> > *xs,
                        const vector<Vector3d> *us,
                        bool rs, 
                        int is, int ie,
                        int dis, int dit,
                        bool dl)
{
  //  glColor4f(

  // set defaults
  if (is == -1)
    is = 0;
  if (ie == -1)
    ie = xs->size()-1;

  assert(is >= 0 && is <= xs->size()-1 && ie >= 0 && ie <= xs->size()-1);
  assert(is <= ie);

  glDisable(GL_LIGHTING);
  //  Viewer::SetColor(rgba[0], rgba[1], rgba[2], rgba[3]);
  glColor4d(rgba[0], rgba[1], rgba[2], rgba[3]);
  
  glLineWidth(lineWidth);
    glBegin(GL_LINE_STRIP);
  for (int i = is; i <= ie; i+=dit) {
    const Matrix3d &g = (*xs)[i].first;
    glVertex3d(g(0,2), g(1,2), 0);
  }
  glEnd();
  glLineWidth(1);
    
  if (rs) {
    for (int i = 0; i < xs->size(); i += dis) {
      if (us && i < us->size())
        Render(&(*xs)[i], &(*us)[i]);
      else
        Render(&(*xs)[i]);
    }
  }
  

  if (dl)
    Render(&xs->back());
  glEnable(GL_LIGHTING);

}
