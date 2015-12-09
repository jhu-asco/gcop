#include "sphereview.h"

using namespace gcop;

SphereView::SphereView(const Sphere &sphere, Matrix4d *g) : 
  Geom3dView("Sphere", g), sphere(sphere)
{
  qobj = gluNewQuadric();
}


SphereView::~SphereView()
{
  free(qobj);
}

void SphereView::RenderGeom()
{  
  glPushMatrix();
  glTranslated(sphere.o[0], sphere.o[1], 0);
  gluSphere(qobj, sphere.r, 10, 10);
  glPopMatrix();
}
