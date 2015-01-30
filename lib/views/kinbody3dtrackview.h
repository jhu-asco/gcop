#ifndef GCOP_KINBODY3DTRACKVIEW_H
#define GCOP_KINBODY3DTRACKVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "kinbody3dtrack.h"
#include "view.h"


namespace gcop {

  /**
   * Rigid body view. Supports rendering either a single state
   * or a whole trajectory of states.
   */
  class Kinbody3dTrackView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    Kinbody3dTrackView(const Kinbody3dTrack &pg);

    virtual ~Kinbody3dTrackView();

    virtual void Render();
    
    virtual bool RenderFrame(int i);
  
    const Kinbody3dTrack &pg;

    float rgba[4];

    bool drawLandmarks;
    bool drawForces;
    double forceScale;


    GLUquadricObj *qobj;   
  };
}

#endif
