#ifndef GCOP_SPHEREVIEW_H
#define GCOP_SPHEREVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "geom3dview.h"
#include "sphere.h"

namespace gcop {

  using namespace Eigen;

  /**
   * Sphere view.
   */
  class SphereView : public Geom3dView {
  public:
    
    SphereView(const Sphere &sphere, Matrix4d *g = 0);
    virtual ~SphereView();

    void RenderGeom();

    const Sphere &sphere;
    GLUquadricObj *qobj;
  };
}

#endif
