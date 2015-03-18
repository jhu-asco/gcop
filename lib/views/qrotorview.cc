#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "qrotorview.h"
#include "viewer.h"
#include "so3.h"
#include "utils.h"

using namespace gcop;
using namespace Eigen;

QrotorView::QrotorView(const Qrotor &sys, 
                       vector< pair<Matrix3d, Vector9d> > *xs,
                       vector< Vector4d > *us) : 
  Body3dView(sys, xs, us), sys(sys)
{
}

void QrotorView::Render(const pair<Matrix3d, Vector9d> *x,
                        const Vector4d *u)
{
  //   glColor4f(1,0.5,0.5,0.5);

  glPushMatrix();
  Transform(x->first, x->second.head<3>());

  //  Viewer::SetColor(0,0.5,0.5,0);

  glColor4dv(rgba);

  glutSolidSphere(.25*sys.l, 10, 10);

  for (int i = 0; i < 4; ++i) {

    if (i==0)
      glColor3f(0,.5,.5);
    else
      glColor3f(0, 0, 0);

    glColor4dv(rgba);
    //    Viewer::SetColor(.5,.5,.5,0);

    glPushMatrix();
    glRotated(90*i, 0, 0, 1);

    glPushMatrix();
    glRotated(90, 0, 1, 0);
    gluCylinder(qobj, sys.l/16, sys.l/16, sys.l, 10, 10);
    glPopMatrix();
    

    glTranslated(sys.l,0,0);
    gluCylinder(qobj, sys.l/8, sys.l/8, sys.l/8, 10, 10);

    //    Viewer::SetColor(0.7,0.7,0.7,0);
    glColor4f(.3,.3,.3, .2);
    glTranslated(0,0,sys.l/8);
    
    //  glRotated(RAD2DEG(ar), 1, 0, 0);
    //  glRotated(RAD2DEG(ap), 0, 1, 0);
    
    //    Viewer::SetColor(0,0,0,0);
    glColor4f(0,0,0,0);
    int i0 = (int)round(rand()/(double)RAND_MAX);
    for (int j=0;j<4;++j)
      gluPartialDisk(qobj, 0, sys.r, 10, 10, i0*45 + 22.5+j*90, 22.5);
    //gluCylinder(tprop, sys.rb, sys.rb, .05, 10, 10);
    
    glPopMatrix();    
  }

  glPopMatrix();
}
