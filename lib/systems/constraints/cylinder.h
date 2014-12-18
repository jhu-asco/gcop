#ifndef GCOP_CYLINDER_H
#define GCOP_CYLINDER_H

#include "constraint.h"

namespace gcop {
  
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic>
    class Cylinder : public Constraint<T, _nx, _nu, _np, 1> {
  public:

  typedef Matrix<double, 1, 1> Vectorgd;
  typedef Matrix<double, 1, _nx> Matrixgxd;
  typedef Matrix<double, 1, _nu> Matrixgud;
  typedef Matrix<double, 1, _np> Matrixgpd;

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;


  Cylinder(const Vector3d &o, double r, double h, double cr = .5);
  
  bool operator()(Vectorgd &g,
                  double t, const T &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrixgxd *dgdx = 0, Matrixgud *dgdu = 0,
                  Matrixgpd *dgdp = 0);    
  
  Vector3d o;      ///< origin of base
  double r;        ///< radius
  double h;        ///< height
  double cr;       ///< collision radius
  
  };
  
  
  template <typename T, int _nx, int _nu, int _np> 
    Cylinder<T, _nx, _nu, _np>::Cylinder(const Vector3d &o, double r, double h, double cr) :
    Constraint<T, _nx, _nu, _np, 1>(), o(o), r(r), h(h), cr(cr)
  {
  }
  
  template <typename T, int _nx, int _nu, int _np> 
    bool Cylinder<T, _nx, _nu, _np>::operator()(Vectorgd &g,
                                             double t, const T &x, const Vectorcd &u,
                                             const Vectormd *rho, 
                                             Matrixgxd *dgdx, Matrixgud *dgdu,
                                             Matrixgpd *dgdp)
    {
      const Vector3d &p = x.second.head(3); // position
      
      Vector2d v = p.head<2>() - o.head<2>();

      double d = v.norm() - cr;  // distance from center of cylinder to system boundary
      
      if (d < 0) {
        cout << "ERR: already colliding" << endl;
      }

      g[0] = r - d; // must be negative for non-collision
      
      if (dgdx) {
        dgdx->Zero();
        if (g[0] > 0)
          dgdx->head(2) = -v/v.norm();
        else
          dgdx->head(2) = v/v.norm();
      }
    }
};



#endif
