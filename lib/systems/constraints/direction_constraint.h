#ifndef GCOP_DIRECTIONCONSTRAINT_H
#define GCOP_DIRECTIONCONSTRAINT_H

#include "constraint.h"
#include "body3dmanifold.h"

namespace gcop {


  /*
    Enforces a constraint between a given body-fixed direction
    vector and the vector between the current position and a given point.
    Currently, it ignores the z-coordinate.
   */
  
template <int nu = 6>
class DirectionConstraint : public Constraint<Body3dState, 12, nu, Dynamic, 3> {
  public:

  using Vectorcd = typename Constraint<Body3dState, 12, nu, Dynamic, 3>::Vectorcd;
  using Matrixgnd = typename Constraint<Body3dState, 12, nu, Dynamic, 3>::Matrixgnd;
  using Matrixgcd = typename Constraint<Body3dState, 12, nu, Dynamic, 3>::Matrixgcd;
  using Matrixgmd = typename Constraint<Body3dState, 12, nu, Dynamic, 3>::Matrixgmd;
  
  DirectionConstraint(const Vector3d& direction,
                      const Vector3d& point);
  
  bool operator()(Vector3d &g,
                  double t, const Body3dState &x, const Vectorcd &u,
                  const VectorXd *p = 0, 
                  Matrixgnd *dgdx = 0, Matrixgcd *dgdu = 0,
                  Matrixgmd *dgdp = 0);
  
  const Vector3d direction;  ///< body-fixed direction
  const Vector3d point;      ///< target point
  
  };
  
    template <int nu>
    DirectionConstraint<nu>::DirectionConstraint(const Vector3d& direction,
                                                     const Vector3d& point) :
    Constraint<Body3dState, 12, nu, Dynamic, 3>(3),//, true),
    direction(direction),
    point(point)
  {
  }
  
    template <int nu>
    bool DirectionConstraint<nu>::operator()(Vector3d &g,
                                                 double t,
                                                 const Body3dState &x,
                                                 const Vectorcd &u,
                                                 const VectorXd *rho, 
                                                 Matrixgnd *dgdx, Matrixgcd *dgdu,
                                                 Matrixgmd *dgdp)
    {
      Vector3d spatial_direction = x.R*direction;
      Vector3d dp = point - x.p;
      double dist = dp.norm();
      if (dist < 1e-3) {
        // TODO(marin): better way to handle singularity near the end?
        dp = spatial_direction;
      } else {
        dp /= dist;
      }
      g = spatial_direction - dp;
      g[2] = 0; // ignore z
      return true;
    }
};


#endif
