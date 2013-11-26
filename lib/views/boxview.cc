#include "boxview.h"

using namespace gcop;

BoxView::BoxView(const Vector3d &d, Matrix4d *g) : 
  Geom3dView("Geom3d View", g), d(d) 
{
}

BoxView::BoxView() : 
  Geom3dView("Geom3d View"), d(1,1,1) 
{
}


void BoxView::RenderGeom()
{  
  glPushMatrix();
  glColor3f(0,1,0);
  glScaled(d[0], d[1], d[2]);
  glutSolidCube(1);
  glPopMatrix();
}
