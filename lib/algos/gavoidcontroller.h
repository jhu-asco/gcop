#ifndef GCOP_GAVOIDCONTROLLER_H
#define GCOP_GAVOIDCONTROLLER_H

#include "controller.h"
#include "body3dmanifold.h"
#include "constraint.h"

namespace gcop {

using namespace std;
using namespace Eigen;

 typedef Constraint<Body3dState, 12, 6, Dynamic, 1> Body3dConstraint;
 typedef Matrix<double, 6, 1> Vector6d;

/**
 * Gyroscopic obstacle avoidance controller, i.e. it produces forces
 * that "rotate" the velocity vector, so that no kinetic energy is inserted into the system
 * Stability is guaranteed for convex obstacles
 *
 *  Author: Marin Kobilarov, 2007 (originally in the DGC library)
 */
class GavoidController : public Controller<Body3dState, Vector6d> {
 public:
  
  /**
   * Gyroscopic forcing controller
   *
   * @param con constraint giving the obstacle distance function
   */
  GavoidController(Body3dConstraint &con);
  
  virtual ~GavoidController();
  
  virtual bool Set(Vector6d &u, double t, const Body3dState &x);
  
  Body3dConstraint &con; ///< obstacle constraint
  
  double sr; ///< sensing radius
  double k;  ///< obstacle avoidance gain
  double kb; ///< breaking gain (optional, zero by default)
  
 };
  
}

#endif
