#ifndef GCOP_DISKCONSTRAINT_H
#define GCOP_DISKCONSTRAINT_H

#include "constraint.h"
#include "disk.h"

namespace gcop {
  
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic>
    class DiskConstraint : public Constraint<T, _nx, _nu, _np, 1> {
  public:

  typedef std::function< int(Vector2d&, const T& ) > ToVector2d;

  typedef Matrix<double, 1, 1> Vectorgd;
  typedef Matrix<double, 1, _nx> Matrixgxd;
  typedef Matrix<double, 1, _nu> Matrixgud;
  typedef Matrix<double, 1, _np> Matrixgpd;

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;


  DiskConstraint(const Disk &disk, double cr = .5);
  
  bool operator()(Vectorgd &g,
                  double t, const T &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrixgxd *dgdx = 0, Matrixgud *dgdu = 0,
                  Matrixgpd *dgdp = 0);
  
  const Disk &disk;   ///< disk
  double cr;          ///< collision radius

  ToVector2d func;

  };
  
  
  template <typename T, int _nx, int _nu, int _np> 
    DiskConstraint<T, _nx, _nu, _np>::DiskConstraint(const Disk& disk, double cr) :
    Constraint<T, _nx, _nu, _np, 1>(), disk(disk), cr(cr)
  {
  }
  
  template <typename T, int _nx, int _nu, int _np> 
    bool DiskConstraint<T, _nx, _nu, _np>::operator()(Vectorgd &g,
                                                      double t, const T &x, const Vectorcd &u,
                                                      const Vectormd *rho, 
                                                      Matrixgxd *dgdx, Matrixgud *dgdu,
                                                      Matrixgpd *dgdp)
    {
      int gi = 0; // index where gradient should start

      Vector2d p;
      if (!std::is_same<T, Vector2d>::value) {
        gi = func(p, x);
      } else {
        p = (Vector2d&)x;
      }

      //      const Vector2d &p = x.head(2); // position
      
      Vector2d v = p - disk.o;

      double d = v.norm() - cr;  // distance from center of disk to system boundary
      
//      if (d < 0) {
//        cout << "[W] Disk:(): in collision, distance is d=" << d << endl;
//      }

      g[0] = disk.r - d; // must be negative for non-collision
      
      if (dgdx) {
        dgdx->Zero();
        if (g[0] > 0)
          dgdx->segment(gi, 2) = -v/v.norm();
        else
          dgdx->segment(gi, 2) = v/v.norm();
      }
    }
};


#endif
