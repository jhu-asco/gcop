#include "cylinderview.h"

using namespace gcop;

CylinderView::CylinderView(double r, double h, Matrix4d *g) : 
  Geom3dView("Cylinder", g), r(r), h(h)
{
  qobj = gluNewQuadric();
}


CylinderView::~CylinderView()
{
  free(qobj);
}

void CylinderView::RenderGeom()
{  
  glPushMatrix();
  glTranslated(0,0,-h/2);
  gluCylinder(qobj, r, r, h, 20, 20);
  glPopMatrix();
}
