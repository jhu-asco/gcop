#ifndef GCOP_SPHERECONSTRAINT_H
#define GCOP_SPHERECONSTRAINT_H

#include "constraint.h"
#include "sphere.h"

namespace gcop {
  
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic>
    class SphereConstraint : public Constraint<T, _nx, _nu, _np, 1> {
  public:

  typedef std::function< int(Vector3d&, const T& ) > ToVector3d;

  typedef Matrix<double, 1, 1> Vectorgd;
  typedef Matrix<double, 1, _nx> Matrixgxd;
  typedef Matrix<double, 1, _nu> Matrixgud;
  typedef Matrix<double, 1, _np> Matrixgpd;

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;


  SphereConstraint(const Sphere &sphere, double cr = .5);
  
  bool operator()(Vectorgd &g,
                  double t, const T &x, const Vectorcd &u,
                  const Vectormd *p = 0, 
                  Matrixgxd *dgdx = 0, Matrixgud *dgdu = 0,
                  Matrixgpd *dgdp = 0);
  
  const Sphere &sphere;   ///< sphere
  double cr;              ///< collision radius
  
  ToVector3d func;        ///< function converting from a generic state T to Vector3d 
  
  };
  
  
  template <typename T, int _nx, int _nu, int _np> 
    SphereConstraint<T, _nx, _nu, _np>::SphereConstraint(const Sphere& sphere, double cr) :
    Constraint<T, _nx, _nu, _np, 1>(), sphere(sphere), cr(cr)
  {
  }
  
  template <typename T, int _nx, int _nu, int _np> 
    bool SphereConstraint<T, _nx, _nu, _np>::operator()(Vectorgd &g,
                                                        double t, const T &x, const Vectorcd &u,
                                                        const Vectormd *rho, 
                                                        Matrixgxd *dgdx, Matrixgud *dgdu,
                                                        Matrixgpd *dgdp)
    {
      int gi = 0; // index where gradient should start

      Vector3d p;
      if (!std::is_same<T, Vector3d>::value) {  
        gi = func(p, x);
      } else {
        p = (Vector3d&)x; // if state is already Vector3d then just copy directly
      }

      Vector3d v = p - sphere.o;

      double d = v.norm() - cr;  // distance from center of sphere to system boundary
      
      if (d < 0) {
        cout << "[W] Sphere:(): in collision, distance is d=" << d << endl;
      }

      g[0] = sphere.r - d; // must be negative for non-collision
      
      if (dgdx) {
        dgdx->Zero();
        if (g[0] > 0)
          dgdx->segment(gi, 3) = -v/v.norm();
        else
          dgdx->segment(gi, 3) = v/v.norm();
      }
    }
};


#endif
