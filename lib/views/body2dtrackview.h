#ifndef GCOP_BODY2DTRACKVIEW_H
#define GCOP_BODY2DTRACKVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "body2dtrack.h"
#include "view.h"


namespace gcop {

  /**
   * Rigid body view. Supports rendering either a single state
   * or a whole trajectory of states.
   */
  class Body2dTrackView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    Body2dTrackView(const Body2dTrack &pg);

    virtual ~Body2dTrackView();

    virtual void Render();
    
    virtual bool RenderFrame(int i);
  
    const Body2dTrack &pg;

    float rgba[4];

    bool drawLandmarks;
    bool drawForces;
    double forceScale;


    GLUquadricObj *qobj;   
  };
}

#endif
