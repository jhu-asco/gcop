#ifndef GCOP_POSEGRAPH2DVIEW_H
#define GCOP_POSEGRAPH2DVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "posegraph2d.h"
#include "view.h"


namespace gcop {

  /**
   * Rigid body view. Supports rendering either a single state
   * or a whole trajectory of states.
   */
  class Posegraph2dView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    Posegraph2dView(const Posegraph2d &pg);

    virtual void Render();
    
    virtual bool RenderFrame(int i);
  
    const Posegraph2d &pg;
    
    float rgba[4];
  };
}

#endif
