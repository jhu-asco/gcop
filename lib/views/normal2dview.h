#ifndef GCOP_NORMAL2DVIEW_H
#define GCOP_NORMAL2DVIEW_H

#include "normal.h"
#include "systemview.h"

namespace gcop {
  using namespace Eigen;

  class Normal2dView : public SystemView<Normal> {
  public:
    /**
     *  Create a normal view of trajectory traj
     * @param sys normal
     * @param xs trajectory
     */
    Normal2dView(vector<Normal> *gds = 0);

    virtual ~Normal2dView();       

    void Render(const Normal &gd);
    
    void Render(const vector<Normal> &gds,
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);

    bool wire;
    bool texture;
    double s; ///< height scaling (default is 1)
  };
}

#endif
