#ifndef GCOP_SHELL_H
#define GCOP_SHELL_H

#include "constraint.h"

namespace gcop {
  
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic>
    class Shell : public Constraint<T, _nx, _nu, _np, 1> {
  public:

  typedef Matrix<double, 1, 1> Vectorgd;
  typedef Matrix<double, 1, _nx> Matrixgxd;
  typedef Matrix<double, 1, _nu> Matrixgud;
  typedef Matrix<double, 1, _np> Matrixgpd;

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;


  Shell(const Vector3d &o, double r, double cr = .5);
  
  bool operator()(Vectorgd &g,
                  double t, const T &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrixgxd *dgdx = 0, Matrixgud *dgdu = 0,
                  Matrixgpd *dgdp = 0);    
  
  Vector3d o;      ///< origin
  double r;        ///< radius
  double cr;       ///< collision radius
  
  };
  
  
  template <typename T, int _nx, int _nu, int _np> 
    Shell<T, _nx, _nu, _np>::Shell(const Vector3d &o, double r, double cr) :
    Constraint<T, _nx, _nu, _np, 1>(), o(o), r(r), cr(cr)
  {
  }
  
  template <typename T, int _nx, int _nu, int _np> 
    bool Shell<T, _nx, _nu, _np>::operator()(Vectorgd &g,
                                             double t, const T &x, const Vectorcd &u,
                                             const Vectormd *rho, 
                                             Matrixgxd *dgdx, Matrixgud *dgdu,
                                             Matrixgpd *dgdp)
    {
      const Vector3d &p = x.second.head(3); // position
      
      Vector3d v = p - o;

      double d = v.norm() - cr;  // distance from center of shell to system boundary
      
      if (d < 0) {
        cout << "ERR: already colliding" << endl;
      }

      g[0] = r - d; // must be negative for non-collision
      
      if (dgdx) {
        dgdx->Zero();
        if (g[0] > 0)
          dgdx->segment(3,3) = -v/v.norm();
        else
          dgdx->segment(3,3) = v/v.norm();
      }
    }
};



#endif
