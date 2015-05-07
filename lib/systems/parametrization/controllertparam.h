// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_CONTROLLERTPARAM_H
#define GCOP_CONTROLLERTPARAM_H

#include "tparam.h"
#include "controller.h"
#include "constraint.h"
#include <assert.h>
#include "utils.h"
#include "sampler.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  /**
   * Trajectory parametrization using a parametrized controller (a parametric policy).
   * The parameters are typically the gains of the controller, but could also be
   * any other policy parameters.
   *
   * Author: Marin Kobilarov (c) 2005-2015
   */
  template <typename Tx, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic,
    int _ntp = Dynamic,
    typename Tc = VectorXd> class ControllerTparam : public Tparam<Tx, nx, nu, np, _ntp> {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;
    
    typedef Matrix<double, np, 1> Vectormd;
    
    typedef Matrix<double, _ntp, 1> Vectorntpd;
    typedef Matrix<double, _ntp, _ntp> Matrixntpd;
    
  public:
    ControllerTparam(System<Tx, nx, nu, np> &sys, 
                     Controller<Tx, Vectorcd, Vectorntpd, Tc> &ctrl,
                     int ntp = _ntp,
                     Constraint<Tx, nx, nu, np, 1> *con = 0);
    
    void To(Vectorntpd &s, 
            const vector<double> &ts, 
            const vector<Tx> &xs, 
            const vector<Vectorcd> &us,
            const Vectormd *p = 0);
    
    void From(vector<double> &ts, 
              vector<Tx> &xs, 
              vector<Vectorcd> &us,
              const Vectorntpd &s,
              Vectormd *p = 0);
    
    Controller<Tx, Vectorcd, Vectorntpd, Tc> &ctrl; ///< controller
    Constraint<Tx, nx, nu, np, 1> *con;         ///< scalar-valued path constraint to check feasibility
 
    bool stoch;                                ///< whether to simulate process noise (true by default)
    Sampler<Tc> *contextSampler;                ///< some stochastic problems might require sampling the context
  };
  
  template <typename Tx, int nx, int nu, int np, int _ntp, typename Tc> 
    ControllerTparam<Tx, nx, nu, np, _ntp, Tc>::ControllerTparam(System<Tx, nx, nu, np> &sys, Controller<Tx, Vectorcd, Vectorntpd, Tc> &ctrl, int ntp, Constraint<Tx, nx, nu, np, 1> *con) :  Tparam<Tx, nx, nu, np, _ntp>(sys, ntp), ctrl(ctrl), con(con), stoch(true), contextSampler(0) {
    assert(ntp > 0);
  }

  template <typename Tx, int nx, int nu, int np, int _ntp, typename Tc> 
    void ControllerTparam<Tx, nx, nu, np, _ntp, Tc>::To(Vectorntpd &s,
                                                        const vector<double> &ts, 
                                                        const vector<Tx> &xs, 
                                                        const vector<Vectorcd> &us,
                                                        const Vectormd *p) {
    cout << "[W] ControllerTparam::To: calling this is useless for this type of paremetrization!" << endl;    
  }  
  
  template <typename Tx, int nx, int nu, int np, int _ntp, typename Tc> 
    void ControllerTparam<Tx, nx, nu, np, _ntp, Tc>::From(vector<double> &ts, 
                                                          vector<Tx> &xs, 
                                                          vector<Vectorcd> &us,
                                                          const Vectorntpd &s,
                                                          Vectormd *p) {

    if (contextSampler) {
      // sample until successful
      do {
        (*contextSampler)(ctrl.c); 
      } while(!ctrl.SetContext(ctrl.c));
    }

    ctrl.SetParams(s);
    Matrix<double, 1, 1> g;

    // sample noise term, assume constant along path
    Vectornd w;
    if (stoch) 
      for (int j=0; j < w.size(); ++j)
        w[j] = random_normal();

    for (int i = 0; i < us.size(); ++i) {
      ctrl.Set(us[i], ts[i], xs[i]);

      if (stoch)
        this->sys.Step(xs[i+1], ts[i], xs[i], us[i], ts[i+1] - ts[i], w, p);      
      else 
        this->sys.Step(xs[i+1], ts[i], xs[i], us[i], ts[i+1] - ts[i], p);
      
      if (con) { // check for feasibility
        (*con)(g, ts[i+1], xs[i+1], us[i], p);
        if (g[0] > 0) {  // if infeasible then set all remaining states to current state
          for (int j = i+2; j < xs.size(); ++j) {
            xs[j] = xs[i+1];         
            us[j-1].setZero();
          }
          break;
        }
      }
    }
  }
}

#endif
