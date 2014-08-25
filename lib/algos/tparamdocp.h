#ifndef GCOP_TPARAMDOCP_H
#define GCOP_TPARAMDOCP_H

#include "docp.h"
#include "tparam.h"

namespace gcop {

  using namespace std;
  using namespace Eigen;
 
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic,
    int nz = Dynamic,
    int _ntp = Dynamic> class TparamDocp : 
    public Docp<T, nx, nu, np, nz>{
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  
    
  public:
    /**
     * Create an optimal control problem using a system, a cost, and 
     * a trajectory given by a sequence of times, states, and controls. 
     * The times ts must be given, the initial state xs[0] must be set, and
     * the controls us will be used as an initial guess for the optimization.
     *
     * After initialization, every call to Iterate() will optimize the 
     * controls us and states xs and modify them accordingly. Problems involving
     * time-optimization will also modify the sequence of times ts.
     * 
     * @param sys system
     * @param cost cost
     * @param tp trajectory paramtrization
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     */
    TparamDocp(System<T, nx, nu, np, nz> &sys, Cost<T, nx, nu> &cost, 
               Tparam<T, nx, nu, np, nz, _ntp> &tp,
               vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us,
               Vectormd *p = 0);
    
    virtual ~TparamDocp();

    /**
     * Perform one DOCP iteration. Internally calls:
     * are updated. 
     */
    void Iterate();

    Tparam<T, nx, nu, np, nz, _ntp> &tp; ///< trajectory parametrization
  };

  
  template <typename T, int nx, int nu, int np, int nz, int _ntp> 
    TparamDocp<T, nx, nu, np, nz, _ntp>::TparamDocp(System<T, nx, nu, np, nz> &sys, 
                                                    LqCost<T, nx, nu> &cost,
                                                    Tparam<T, nx, nu, np, nz, _ntp> &tp,
                                                    vector<double> &ts, 
                                                    vector<T> &xs, 
                                                    vector<Matrix<double, nu, 1> > &us,
                                                    Matrix<double, np, 1> *p) : 
    Docp<T, nx, nu, np, nz>(sys, cost, ts, xs, us, p, false) {
  }
  
  template <typename T, int nx, int nu, int np, int nz, int _ntp> 
    TparamDocp<T, nx, nu, np, nz, _ntp>::~TparamDocp()
    {
    }  
  
  template <typename T, int nx, int nu> 
    void Tparamdocp<T, nx, nu>::Iterate() {
    
  }
}

#endif
