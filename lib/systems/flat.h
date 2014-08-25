#ifndef GCOP_FLAT_H
#define GCOP_FLAT_H

#include "system.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  /**
   * General flat output definition
   *
   * Author: Marin Kobilarov (c) 2005--2013
   */
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic,
    int nz = Dynamic,
    int _ny = Dynamic> class Flat {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  

    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, _ny, 1> Vectoryd;

  public:
    /**
     * Define flat outputs for system s with state order ox and control order ou
     * @param sys system
     * @param ny flat output dimension
     * @param ox state order
     * @param ou control order
     */
    Flat(System<T, nx, nu, np, nz, _ny> &sys, int ny = _ny, int ox = 1, int ou = 2);
    
    /**
     * Construct state from flat outputs
     * @param x state
     * @param ys a vector of y and its derivatives (y,dy,d2y,...)
     * @return true if paramaters are feasible
     */
    virtual bool State(T& x, 
                       const vector<Vectoryd> &ys);
    
    /**
     * Construct state from flat outputs
     * @param u control
     * @param ys a vector of y and its derivatives (y,dy,d2y,...)
     * @return true if paramaters are feasible
     */
    virtual bool Control(Vectorcd &u, 
                         const vector<Vectormd> &ys);

    System<T, nx, nu, np, nz> &sys;     ///< system

    int ny;                            ///< flat output dimension
    int ox;                            ///< state order
    int ou;                            ///< control order
  };
  
  template <typename T, int nx, int nu, int np, int nz, int _ny> 
    Flat<T, nx, nu, np, nz, _ny>::Flat(System<T, nx, nu, np, nz> &sys, int ny, int ox, int ou) : 
    sys(sys), ny(_ny != Dynamic ? _ny : ny), ox(ox), ou(ou) {
    assert(ny > 0);
    assert(ox > 0);
    assert(ou > 0);
  }

  template <typename T, int nx, int nu, int np, int nz, int _ny> 
    bool Flat<T, nx, nu, np, nz, _ny>::State(T &x,
                                             const vector<Vectoryd> &ys) {
    cout << "[W] Flat::State: unimplemented! Subclasses should override." <<endl;
    return false;
  } 

  template <typename T, int nx, int nu, int np, int nz, int _ny> 
    bool Flat<T, nx, nu, np, nz, _ny>::Control(Vectorcd &u,
                                               const vector<Vectoryd> &ys) {
    cout << "[W] Flat::Control: unimplemented! Subclasses should override." <<endl;
    return false;
  } 
}

#endif
