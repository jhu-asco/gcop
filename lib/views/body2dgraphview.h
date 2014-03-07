#ifndef GCOP_BODY2DGRAPHVIEW_H
#define GCOP_BODY2DGRAPHVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "body2dgraph.h"
#include "view.h"


namespace gcop {

  /**
   * Rigid body view. Supports rendering either a single state
   * or a whole trajectory of states.
   */
  class Body2dGraphView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    Body2dGraphView(const Body2dGraph &pg);

    virtual void Render();
    
    virtual bool RenderFrame(int i);
  
    const Body2dGraph &pg;

    float rgba[4];

    bool drawLandmarks;
   
  };
}

#endif
