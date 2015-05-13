#ifndef GCOP_TPARAM_H
#define GCOP_TPARAM_H

#include "system.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  /**
   * General trajectory parametrization functionality. A trajectory can be parametrizes using
   * a set of discrete controls, a continuous control paramertization, a set of discrete states, 
   * a set of discrete flat outputs, etc...
   *
   * The default implementation is the simplest: a set of discrete controls
   *
   * Author: Marin Kobilarov (c) 2005--2013
   */
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic,
    int _ntp = Dynamic,
    typename Tc = T> class Tparam {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  

    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, _ntp, 1> Vectorntpd;
    
  public:
    Tparam(System<T, nx, nu, np> &sys, int ntp = 0);
    
    /**
     * Convert from paramer s to trajectory (ts,xs,us,p)
     * @param s trajectory parametrization vector
     * @param ts times
     * @param xs states
     * @param us controls
     * @param p system parameters (optional)
     * @return true if conversion was successful
     */
    virtual bool To(Vectorntpd &s, 
                    const vector<double> &ts, 
                    const vector<T> &xs, 
                    const vector<Vectorcd> &us,
                    const Vectormd *p = 0);
    
    /**
     * Convert from trajectory (ts,xs,us,p) to parameters s
     * @param ts times
     * @param xs states
     * @param us controls
     * @param s trajectory parametrization vector
     * @param p system parameters (optional)
     * @return true if conversion was successful
     */
    virtual bool From(vector<double> &ts, 
                      vector<T> &xs, 
                      vector<Vectorcd> &us,
                      const Vectorntpd &s,
                      Vectormd *p = 0);

    virtual bool SetContext(const Tc &c) { return true; }
    
    System<T, nx, nu, np> &sys;     ///< system
    int ntp;                        ///< parameter vector dimension
  };
  
  template <typename T, int nx, int nu, int np, int _ntp, typename Tc> 
    Tparam<T, nx, nu, np, _ntp, Tc>::Tparam(System<T, nx, nu, np> &sys, int ntp) : 
    sys(sys), ntp(_ntp != Dynamic ? _ntp : ntp) {
    assert(ntp > 0);
    std::cout<<"ntp: "<<ntp<<std::endl;//#DEBUG
  }

  template <typename T, int nx, int nu, int np, int _ntp, typename Tc> 
    bool Tparam<T, nx, nu, np, _ntp, Tc>::To(Vectorntpd &s, 
                                         const vector<double> &ts, 
                                         const vector<T> &xs, 
                                         const vector<Vectorcd> &us,
                                         const Vectormd *p) {
    assert(ntp == us.size()*sys.U.n);
    for (int i = 0; i < us.size(); ++i) {
      memcpy(s.data() + i*sys.U.n, us[i].data(), sys.U.n*sizeof(double));
    }
    return true;
    // assume parameter p does not play a role
  }
  
  template <typename T, int nx, int nu, int np, int _ntp, typename Tc> 
    bool Tparam<T, nx, nu, np, _ntp, Tc>::From(vector<double> &ts, 
                                           vector<T> &xs, 
                                           vector<Vectorcd> &us,
                                           const Vectorntpd &s,
                                           Vectormd *p) {
    assert(ntp == us.size()*sys.U.n);
    sys.Reset(xs[0],ts[0]);
    for (int i = 0; i < us.size(); ++i) {
      memcpy(us[i].data(), s.data() + i*sys.U.n, sys.U.n*sizeof(double));
      sys.Step(xs[i+1], us[i], ts[i+1] - ts[i], p);
    }
    return true;
  }
}

#endif
