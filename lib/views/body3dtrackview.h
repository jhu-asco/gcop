#ifndef GCOP_BODY3DTRACKVIEW_H
#define GCOP_BODY3DTRACKVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "body3dtrack.h"
#include "view.h"


namespace gcop {

  /**
   * Rigid body view. Supports rendering either a single state
   * or a whole trajectory of states.
   */
  class Body3dTrackView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    Body3dTrackView(const Body3dTrack &pg);

    virtual ~Body3dTrackView();

    virtual void Render();
    
    virtual bool RenderFrame(int i);
  
    const Body3dTrack &pg;

    float rgba[4];

    bool drawLandmarks;
    bool drawForces;
    double forceScale;


    GLUquadricObj *qobj;   
  };
}

#endif
