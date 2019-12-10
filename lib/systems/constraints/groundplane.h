#ifndef GCOP_GROUNDPLANE_H
#define GCOP_GROUNDPLANE_H

#include "constraint.h"

//A cost for use with GCOP.  Implements a planar constraint above a certain height.
//Ensures that the system "T" has a height above h+cr, 
//where h is the plane height, and cr is the collision radius

namespace gcop {
  
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic>
    class GroundPlane : public Constraint<T, _nx, _nu, _np> {
  public:

  typedef Matrix<double, Dynamic, 1> Vectorgd;
  typedef Matrix<double, Dynamic, _nx> Matrixgxd;
  typedef Matrix<double, Dynamic, _nu> Matrixgud;
  typedef Matrix<double, Dynamic, _np> Matrixgpd;

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;


  GroundPlane(double h, double cr = .5);
  
  bool operator()(Vectorgd &g,
                  double t, const T &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrixgxd *dgdx = 0, Matrixgud *dgdu = 0,
                  Matrixgpd *dgdp = 0);    
  
  double h;        ///< height
  double cr;       ///< collision radius
  
  };
  
  
  template <typename T, int _nx, int _nu, int _np> 
    GroundPlane<T, _nx, _nu, _np>::GroundPlane(double h, double cr) :
    Constraint<T, _nx, _nu, _np>(1), h(h), cr(cr)
  {
    this->ng = 1;
  }
  
  template <typename T, int _nx, int _nu, int _np> 
    bool GroundPlane<T, _nx, _nu, _np>::operator()(Vectorgd &g,
                                             double t, const T &x, const Vectorcd &u,
                                             const Vectormd *rho, 
                                             Matrixgxd *dgdx, Matrixgud *dgdu,
                                             Matrixgpd *dgdp)
    {
      g.resize(1);
      const Vector3d &p = x.p; // position
      
      double v = p(2) - h;

      double d = v - cr;  // distance from center of cylinder to system boundary
      
      g[0] = -d; // must be negative for non-collision
      
      if (dgdx) {
        dgdx->resize(1, _nx);
        dgdx->setZero();
        //THIS IS A BUG.  For states that are not Body3dState, it won't necessarily work.
        //It assumes that the position indices of the x vector are 3,4,5, but other states use 0,1,2
        //DO NOT TRUST 
        if (std::is_same<T, Body3dState>::value){
          dgdx->row(0)[5] = -1;
        } else {
          dgdx->row(0)[2] = -1;
        }
      }
    }
};



#endif
