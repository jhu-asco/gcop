#ifndef GCOP_CYLINDERVIEW_H
#define GCOP_CYLINDERVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "geom3dview.h"


namespace gcop {

  using namespace Eigen;

  /**
   * Cylinder view.
   */
  class CylinderView : public Geom3dView {
  public:
    
    CylinderView(double r = .5, double h = 1, Matrix4d *g = 0);
    virtual ~CylinderView();

    void RenderGeom();

    double r;
    double h;
    GLUquadricObj *qobj;
  };
}

#endif
