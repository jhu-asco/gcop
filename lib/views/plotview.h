#ifndef GCOP_PLOTVIEW_H
#define GCOP_PLOTVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "view.h"
#include <Eigen/Dense>


namespace gcop {

  using namespace Eigen;

  /**
   * Discrete Dynamical plot view.
   */
  class PlotView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    PlotView(const char *name = 0, Matrix4d *g = 0);

    virtual ~PlotView();

    virtual void Render();
    
    

    /**
     * Set the rgb color of the trajectory
     * @param rgb rgb vector
     */
    void SetColor(const double rgba[4]);

    double rgba[4];      ///< color
  };
}

#endif
