#ifndef GCOP_POINT3DCONTROLLER_H
#define GCOP_POINT3DCONTROLLER_H

#include "controller.h"
#include "point3d.h"

namespace gcop {

  using namespace std;
  using namespace Eigen;

  /**
   * Basic PD controller for a point-mass system in 3d
   *
   * Author: Marin Kobilarov
   */
  class Point3dController : public Controller< Point3dState, Vector3d> {
  public:

    /**
     * Point3D PD controller
     *
     * @param sys multi-body system
     * @param xd desired state (optional, set to origin by default)
     * @param ad desired acceleration (optional, set to zero by default)
     */
    Point3dController(const Point3d &sys,
                      Point3dState *xd = 0, 
                      Vector3d *ad = 0);
    
    virtual void Set(Vector3d &u, double t, const Point3dState &x);

    virtual ~Point3dController();

    const Point3d &sys; ///< point mass system

    Point3dState *xd;   ///< desired state (origin by default)
    Vector3d *ad;       ///< desired acceleration (zero by default)

    Vector3d Kp;    ///< proportional terms  (ones by default)
    Vector3d Kd;    ///< derivative  terms  (ones by default)    
  };
};

#endif
