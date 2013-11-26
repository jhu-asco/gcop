#ifndef GCOP_PARTICLE2DVIEW_H
#define GCOP_PARTICLE2DVIEW_H

#include "particle2d.h"
#include "systemview.h"

namespace gcop {
  using namespace Eigen;

  class Particle2dView : public SystemView<Vector4d> {
  public:
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    Particle2dView(const Particle2d &sys,
                   vector<Vector4d> *xs = 0);

    virtual ~Particle2dView();
       

    void Render(const Vector4d &x);
    
    void Render(const vector<Vector4d> &xs, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Particle2d &sys;
  };
}

#endif
