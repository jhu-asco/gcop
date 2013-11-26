#ifndef GCOP_HELI_H
#define GCOP_HELI_H

#include "body3d.h"

namespace gcop {
  
  /**
   * Idealized small helicopter modelled as a rigid body
   * with 6 degrees of freedom and 6 velocities and 4 control variables.
   * The heli position is denoted \f$ p \in \mathbb{R}^3 \f$, orientation
   * \f$ R \in SO(3)\f$ parametrized using roll \f$ \phi \f$, pitch \f$ \theta \f$, 
   * and yaw \f$ \psi \f$. The linear velocity is denoted \f$ v=(v_x,v_y,v_z) \in \mathbb{R}^3 \f$, 
   * and the angular velocity \f$ \omega=(\omega_x,\omega_y, \omega_z) \in \mathbb{R}^3 \f$. The whole state of the helicopter is denoted by the 12-tuple 
   * \f$ \mathbf{x} = (\omega_x,\omega_y, \omega_z, v_x, v_y, v_z, \phi, \theta, \psi, p_x, p_y, p_z) \f$. 

   * The controls consist of the angles that control the main rotor by pitching it
   * forward (pitch \f$ \gamma_p \f$) and tilting it sideways (roll \f$ \gamma_r \f$)
   * 2 remaining inputs are the collective and (yaw)rudder forces
   * so that \f$ \mathbf{u} = (u_c, u_y, \gamma_p, \gamma_r) \f$. 
   *
   * It is assumed that the center of mass (COM) lies directly below the main rotor, i.e. that
   * the main rotor axis is aligned with the body fixed z-axis (with origin at the COM).
   * Similar assumption is made about the rear rotor application point with respect to the
   * negative body-fixed x-axis. The rear rotor produces rudder force in the y-axis direction.
   * The default heli parameters are distance b/n COM and top rotor\f$ r_t = .3 \f$
   * and distance b/n COM and rear rotor \f$ r_b = 1 \f$, mass \f$ 5 kg.\f$.  Moments
   * of inertia can be taken from a comparable small RC helicopter but here we assume
   * that the helicopter mass is concentrated in a sphere of radius rt/2 and mass mb fixed
   * at the COM. Rotational moments of inertia are computed using the sphere as well as a 
   * a virtual disk or radius rb and weight mr 
   * representing the spinning top rotor fixed at distance rt above the COM
   * perpendicular to the body-fixed z-axis.
   *
   * 
   * Denote the whole configuration by \f$ g \in SE(3) \f$, the whole velocity by
   * \f$ \xi \in \mathfrak{se}(3)\f$, defined by
   * \f{align*} 
   * g = \left[\begin{array}{cc} 
   * R & p \\ 
   * \mathbf{0} & 1 
   * \end{array}\right], \quad \xi = \left[ \begin{array}{cc}
   * \widehat\omega & v \\
   * \mathbf{0} & 0 
   * \end{array}\right].
   * \f}
   * 
   * Denote the transformation corresponding to
   * roll and pitch only by \f$ g(\phi,\theta) \in SE(3) \f$.
   * Denote a rotation corresponding to roll \f$\phi\f$, and pitch \f$\theta\f$ by
   * \f$ R(\phi, \theta) \in SO(3) \f$. 
   *
   * Denoting the vectors 
   * \f$ f_c = R(\gamma_r, \gamma_p) (0,0,u_c)\in \mathbb{R}^3 \f$ and
   * \f$ f_y = (0, -u_y, 0) \in \mathbb{R}^3\f$.
   * 
   * The control forces \f$ f_u \in \mathfrak{se}(3)^*\f$ 
   * acting on the system in its body frame are 
   *
   * \f{align*} 
   * f_u = ((0,0,r_t) \times f_c, f_c) + ((-r_b,0,0) \times f_y, f_y).
   * \f}
   *
   * Now let's check what conditions are necessary for motion invariant
   * to translations and rotations around the \f$ z \f$-axis, i.e. to 
   * \f$ G^{\prime} = SE(2)\times \mathbb{R}\f$ transformations. 
   * This means that there is a velocity 
   * \f$ \xi^{\prime}\in \mathfrak{se}(2) \times \mathbb{R} \f$, that can 
   * be expressed as \f$ \xi^{\prime}=(0, 0, \omega_z, v_x, v_y, v_z)\f$, for
   * which \f$ \xi = \operatorname{Ad}_{g(\phi,\theta)^{-1}}\xi^{\prime} \f$ 
   * is a relative equilibria for the whole system on \f$ SE(3) \f$, i.e. 
   * \f$ \dot \xi = 0 \f$ and \f$ g(t) = g(0)\exp(t\xi) \f$.
   * 
   * This velocity is obtained by satisfying the invariance conditions
   *
   * \f{align*} 
   * \operatorname{ad}_\xi^* \operatorname{\mathbb{I}} \xi + f_u + \operatorname{Ad}^{\ast}_{g(\phi,\theta)} f_{ext} = 0,
   * \f}
   * using the \f$G^\prime\f$-invariant external force 
   * \f$ f_{ext} = (0,0,0,0,0,f_g m) \in \mathfrak{se}(3)^* \f$ with  
   * the scalar \f$ f_g \f$ denoting acceleration due to gravity 
   * (e.g. on Earth \f$ f_g=-9.81 \f$).
   *
   * These conditions can be simplified if one assumes that the moments 
   * of rotational inertia around the \f$ y \f$ and \f$ z \f$ axis are identical
   * In this case the invariance requires that:
   *
   * \f{align*}
   * & \theta = 0, \quad u_y = 0, \quad \gamma_p = 0, \quad \gamma_r = 0, \\
   * & \phi = \arctan(-w_z v_x/f_g),\\
   * & u_c = m(\cos\phi f_g - w_z v_x \sin\phi).
   * \f}
   * 
   * Author: Marin Kobilarov -- Copyright (C) 2006
   */
  class Heli : public Body3d<4> {
  public:
    Heli();    
    
    double rt;  ///< distance from center of mass to main rotor (top)
    double rb;  ///< distance from center of mass to rear rotor (back)
  };
}

#endif
