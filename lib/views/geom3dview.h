#ifndef GCOP_GEOM3DVIEW_H
#define GCOP_GEOM3DVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "view.h"
#include <Eigen/Dense>


namespace gcop {

  using namespace Eigen;

  /**
   * Discrete Dynamical geom3d view.
   */
  class Geom3dView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    Geom3dView(const char *name = 0, Matrix4d *g = 0);

    virtual ~Geom3dView();

    virtual void RenderGeom() = 0;

    virtual void Render();

    static void Transform(const Matrix4d &g);
    
    /**
     * Set the rgb color of the trajectory
     * @param rgb rgb vector
     */
    void SetColor(const double rgba[4]);

    Matrix4d *g;         ///< global pose (optional)

    Matrix4d g0;         ///< origin (set to Id by default)

    double rgba[4];      ///< color
  };
}

#endif
