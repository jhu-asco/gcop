#ifndef GCOP_LSSENSORCOST_H
#define GCOP_LSSENSORCOST_H

#include <Eigen/Dense>
#include <iostream>
#include "system.h"
#include "manifold.h"
#include "sensorcost.h"

namespace gcop {

  using namespace Eigen;
  using namespace std;
  
  /**
   * Cost interface for optimal control on sensor manifold. 
   * Defines a cost function and means to compute difference between
   * two states on a sensor manifold. 
   *
   * Subclasses should provide implementation of either a regular const function L
   * or a parameter-dependent cost function Lp 
   * (e.g. for sys id / adaptive control / optimal design problems)
   *
   * This defines cost only for a single sensor. Multiple sensors will and multiple systems based cost functions should be adapted somehow
   *
   * Author1: Marin Kobilarov marin(at)jhu.edu
   *
   * Author2: Gowtham Garimella ggarime1(at)jhu.edu
   *
   */
  template <typename T = VectorXd, 
    int _nx = Dynamic,
    int _nu = Dynamic,
    int _np = Dynamic,
    int _ng = Dynamic,
    typename Tz = VectorXd,
    int _nz = Dynamic> class LsSensorCost : public SensorCost<T, _nx, _nu, _np, Tz, _nz> {
    
  public:
  
  typedef Matrix<double, _ng, 1> Vectorgd;
  typedef Matrix<double, _ng, _nx> Matrixgxd;
  typedef Matrix<double, _ng, _nu> Matrixgud;
  typedef Matrix<double, _ng, _np> Matrixgpd;
  typedef Matrix<double, _ng, _nz> Matrixgzd;
  
  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;
  
  typedef Matrix<double, _nx, _nx> Matrixnd;
  typedef Matrix<double, _nx, _nu> Matrixncd;
  typedef Matrix<double, _nu, _nx> Matrixcnd;
  typedef Matrix<double, _nu, _nu> Matrixcd;
  
  //  typedef Matrix<double, Dynamic, 1> Vectormd;
  typedef Matrix<double, _np, _np> Matrixmd;
  typedef Matrix<double, _nx, _np> Matrixnmd;
  typedef Matrix<double, _np, _nx> Matrixmnd;

  typedef Matrix<double, _nz, 1> Vectorrd;
  typedef Matrix<double, _nz, _nz> Matrixrd;
  typedef Matrix<double, _nz, _nx> Matrixrnd;
  typedef Matrix<double, _nz, _nu> Matrixrcd;
  typedef Matrix<double, _nz, _np> Matrixrmd;

  /**
   * create a cost interface for optimal estimation problem
   * @param sys provides the state manifold which is used to perform addition/subtraction of states
   * @param sensor provides the sensor manifold and measurement model for computing the cost
   * @param ng number of residuals per time-step
   */
  LsSensorCost(System<T, _nx, _nu, _np> &sys, Manifold<Tz, _nz> &Z, int ng = 0);

  /**
   * Residual cost with noise as input (For estimation problems) #TODO Modify this to directly accept z sensor estimation instead of x, u
   * @param g residual
   * @param t time
   * @param z sensor measurement
   * @param w noise
   * @param h time-step
   * @param p parameter (optional)
   * @param dgdw derivative wrt w (optional)
   * @param dgdp derivative wrt p (optional)
   * @param dgdz derivative wrt p (optional)
   */
  virtual bool Res(Vectorgd &g, 
                     double t, const Tz &z,
                     const Vectornd &w, const Vectormd &p, double h,
                     Matrixgxd *dgdw = 0, Matrixgpd *dgdp = 0,
                     Matrixgzd *dgdz = 0);    

  /**
   * Residual Cost independent of time and only dependent on parameters
   * @param g residual for parameters
   * @param p parameter 
   * @param dgdp derivative wrt p
   *
   */
  virtual bool Resp(Vectormd &gp, 
                   const Vectormd &p,
                   Matrixmd *dgdp = 0);    

  int ng;       ///< number of residuals per time-step
  Vectorgd g;   ///< residuals per time-step
  Vectormd gp;  ///< residual dependent on parameters
};

  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz> 
    LsSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::LsSensorCost(System<T, _nx, _nu, _np> &sys, Manifold<Tz, _nz> &Z,  int ng) : 
    SensorCost<T, _nx, _nu, _np, Tz, _nz>(sys, Z), ng(_ng != Dynamic ? _ng : ng)
    {
      assert(ng > 0);
      if (_ng == Dynamic)
        g.resize(ng);
      const int np = sys.P.n;
      if(_np == Dynamic)
        gp.resize(np);
    } 

  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz> 
    bool LsSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::Res(Vectorgd &g, 
                     double t, const Tz &z, 
                     const Vectornd &w, const Vectormd &p, double h,
                     Matrixgxd *dgdw, Matrixgpd *dgdp,
                     Matrixgzd *dgdz) {
    cout << "[W] LsSensorCost:Res: unimplemented! Subclasses should override." << endl;
    return false;
  }

  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz> 
    bool LsSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::Resp(Vectormd &g, 
                                             const Vectormd &p,
                                             Matrixmd *dgdp)
    {
      cout << "[W] LsSensorCost:Resp: unimplemented! Subclasses should override." << endl;
      return false;
    }
}

#endif

