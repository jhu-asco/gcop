#include "hrotorgeom3dview.h"
#include "viewer.h"

using namespace gcop;

HrotorGeom3dView::HrotorGeom3dView(const Hrotor &sys,
                                   Matrix4d *g) : 
Geom3dView("Hrotor", g), sys(sys) 
{
  qobj = gluNewQuadric();
}


HrotorGeom3dView::~HrotorGeom3dView()
{
  free(qobj);
}



void HrotorGeom3dView::RenderGeom() 
{
  //  Viewer::SetColor(0,0.5,0.5,0);
  //  glColor3f(0,.5,.5);

  glPushMatrix();
  glTranslated(0,0,-.4*sys.l);
  glutSolidCube(.4*sys.l);
  glPopMatrix();

  glutSolidSphere(.25*sys.l, 10, 10);

    
  //  glDisable(GL_LIGHTING);

  //  glPushMatrix();
  //  glScalef(.5*sys.l,.5*sys.l,.5*sys.l);
  //  Viewer::DrawFrame(this->qobj);
  //  glPopMatrix();

  //  glEnable(GL_LIGHTING);

  for (int i = 0; i < 3; ++i) {    

    if (i==0 || i==2)
      glColor3f(0,.5,.5);
    else
      glColor3f(0, 0, 0);
    
    Viewer::SetColor(1,1,1,0);

    glPushMatrix();
    glRotated(60 + 120*i, 0, 0, 1);

    glPushMatrix();
    glRotated(90, 0, 1, 0);
    gluCylinder(qobj, sys.l/16, sys.l/16, sys.l, 10, 10);
    glPopMatrix();    

    glTranslated(sys.l, 0, -sys.l/8);
    gluCylinder(qobj, sys.l/8, sys.l/8, sys.l/4, 10, 10);

    glColor4f(.3,.3,.3, .2);
    glTranslated(0,0,sys.l/4);
    
    Viewer::SetColor(0,0,0,0);

    int i0 = (int)round(rand()/(double)RAND_MAX);
    for (int j=0;j<4;++j)
      gluPartialDisk(qobj, 0, sys.r, 10, 10, i0*45 + 22.5+j*90, 22.5);

    glTranslated(0,0,-sys.l/4);
    
    i0 = (int)round(rand()/(double)RAND_MAX);
    for (int j=0;j<4;++j)
      gluPartialDisk(qobj, 0, sys.r, 10, 10, i0*45 + 22.5+j*90, 22.5);
    
    glPopMatrix();    
  }
}
