#ifndef GCOP_DYNVISINSVIEW_H
#define GCOP_DYNVISINSVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include <stdlib.h>
#include <string.h>
#include "dynvisins.h"
#include "body3dview.h"


namespace gcop {

  /**
   * Discrete Dynamical visins view.
   */
  class DynVisInsView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    DynVisInsView(const DynVisIns &vi);

    virtual ~DynVisInsView();

    virtual void Render();
    
    virtual bool RenderFrame(int i);
    
    const DynVisIns &vi;

    Body3d<> body;
    Body3dView<> bodyView;
  };
}

#endif
