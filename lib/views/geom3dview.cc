#include "geom3dview.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "viewer.h"
#include "so3.h"
#include "utils.h"

using namespace gcop;

Geom3dView::Geom3dView(const char *name, Matrix4d *g) : 
  View(name ? name : "Geom3d View"), g(g), g0(Matrix4d::Identity()) {
  rgba[0] = 1;
  rgba[1] = 1;
  rgba[2] = 1;
  rgba[3] = 0.5;
}

Geom3dView::~Geom3dView() {
}

void Geom3dView::SetColor(const double rgba[4]) {
  memcpy(this->rgba, rgba, 4*sizeof(double));
}


void Geom3dView::Render()
{
  //  glColor4dv(rgba);
  Viewer::SetColor(rgba[0], rgba[1], rgba[2], 0);

  glPushMatrix();
  if (g)
    Transform((*g)*g0);
  else
    Transform(g0);
    
  RenderGeom();
  glPopMatrix();
}

void Geom3dView::Transform(const Matrix4d &g) 
{  
  glTranslated(g(0,3), g(1,3), g(2, 3));
  
  Vector3d e;
  SO3::Instance().log(e, g.topLeftCorner<3,3>());
  double n = e.norm();
  if (n > SO3::Instance().tol) {
    e = e/n;
    glRotated(RAD2DEG(n), e[0], e[1], e[2]);
  }
}
