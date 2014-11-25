#ifndef GCOP_SENSORCOST_H
#define GCOP_SENSORCOST_H

#include <Eigen/Dense>
#include <iostream>
#include "system.h"
#include "manifold.h"

namespace gcop {

  using namespace Eigen;
  using namespace std;
  
  /**
   * Cost interface for optimal control on sensor manifold. (Estimation Problems)
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
   */
   //T1 can be the input manifold of sensor which need not be the same as the input of system since we are sharing the same sensor between multiple systems
   //typename T1 = T 

  template <typename T = VectorXd, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic,
    typename Tz = VectorXd,
    int _nz = Dynamic> class SensorCost {
  public:
  
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
   * create a cost interface
   * @param sys provides the system manifold X for performing addition and substraction of states
   * @param Z provides sensor manifold for performing addition and substraction of measurements
   */
  SensorCost(System<T, _nx, _nu, _np> &sys, Manifold<Tz, _nz>  &Z);

  /**
   * Time dependent Cost function L_i with noise as a part of cost function (Optimal Estimation)
   * @param t time
   * @param z sensor measurement
   * @param w noise
   * @param h time-step
   * @param Lw derivative wrt w
   * @param Lp derivative wrt p
   * @param Lz derivative wrt z
   * @param Lpp derivative wrt p,p
   * @param Lww derivative wrt w,w
   * @param Lzz derivative wrt z,z
   * @param Lpw derivative wrt p,w
   * @param Lzw derivative wrt z,w
   * @param Lzp derivative wrt z,p
   */
  virtual double L(double t, const Tz &z, const Vectornd &w,
                   const Vectormd &p, double h, 
                   Vectornd *Lw = 0, Matrixnd* Lww = 0,
                   Vectormd *Lp = 0, Matrixmd* Lpp = 0,
                   Vectorrd *Lz = 0, Matrixrd* Lzz = 0,
                   Matrixmnd *Lpw = 0, Matrixrnd* Lzw = 0, Matrixrmd* Lzp = 0);  
  /**
   * Part of cost function which is independent of time and
   * only dependent on parameters
   * The net cost for a trajectory is sum i:1->N [L_i] + Lp
   * @param p parameters
   * @param Lp derivative wrt p
   * @param Lpp derivative wrt p,p
   */
  virtual double Lp(const Vectormd &p,
                   Vectormd *Lp = 0, Matrixmd *Lpp = 0);

  System<T, _nx, _nu, _np> &sys;  ///< system  
  Manifold<Tz, _nz> &Z;  ///< system  
};


  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
    SensorCost<T, _nx, _nu, _np, Tz, _nz>::SensorCost(System<T, _nx, _nu, _np> &sys, 
                                 Manifold<Tz, _nz> &Z) : sys(sys), 
                                                         Z(Z) {
  }  

  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
  double SensorCost<T, _nx, _nu, _np, Tz, _nz>::L(double t, const Tz &z, 
                                             const Vectornd &w,
                                             const Vectormd &p, double h,
                                             Vectornd *Lw, Matrixnd* Lww,
                                             Vectormd *Lp, Matrixmd* Lpp,
                                             Vectorrd *Lz, Matrixrd* Lzz,
                                             Matrixmnd *Lpw, Matrixrnd* Lzw, Matrixrmd* Lzp)
  {
    cout << "[W] Cost:L: unimplemented!" << endl;
    return 0;
  }

  template <typename T, int _nx, int _nu, int _np, typename Tz, int _nz> 
  double SensorCost<T, _nx, _nu, _np, Tz, _nz>::Lp(const Vectormd &p,
                   Vectormd *Lp, Matrixmd *Lpp)
  {
    return 0;
  }
}

#endif

