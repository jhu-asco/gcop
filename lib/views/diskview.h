#ifndef GCOP_DISKVIEW_H
#define GCOP_DISKVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "geom3dview.h"
#include "disk.h"

namespace gcop {

  using namespace Eigen;

  /**
   * Disk view.
   */
  class DiskView : public Geom3dView {
  public:
    
    DiskView(const Disk &disk, Matrix4d *g = 0);
    virtual ~DiskView();

    void RenderGeom();

    const Disk &disk;
    GLUquadricObj *qobj;
  };
}

#endif
