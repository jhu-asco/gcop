#ifndef GCOP_BOXVIEW_H
#define GCOP_BOXVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "geom3dview.h"
#include <Eigen/Dense>


namespace gcop {

  using namespace Eigen;

  /**
   * Box view.
   */
  class BoxView : public Geom3dView {
  public:
    
    BoxView();
    BoxView(const Vector3d &d, Matrix4d *g = 0);

    void RenderGeom();

    Vector3d d;  ///< dimensions
  };
}

#endif
