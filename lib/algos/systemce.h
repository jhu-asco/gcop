// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_SYSTEMCE_H
#define GCOP_SYSTEMCE_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include "system.h"
#include "cost.h"
#include "ce.h"
#include <cmath>
#include "rn.h"
#include <limits>

#include "tparam.h"
#include "creator.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
  /**
   * The CE method for dynamical systems
   * 
   *  Note: the vector dus should denote the variances along the diagonal of the initial distribution
   *
   * Authors: Marin Kobilarov, Gowtham Garimella
   */
  template <typename T, 
    int n = Dynamic, 
    int c = Dynamic, 
    int np = Dynamic, 
    int ntp = Dynamic,
    typename Tc = T> 
    class SystemCe {
    
    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, n, n> Matrixnd;
    typedef Matrix<double, n, c> Matrixncd;
    typedef Matrix<double, c, n> Matrixcnd;
    typedef Matrix<double, c, c> Matrixcd;  

    typedef Matrix<double, ntp, 1> Vectortpd;  
    typedef Matrix<double, ntp, ntp> Matrixtpd;  
    
    //External render Function for trajectories:
    typedef void(RenderFunc)(int, vector<T>&);
  public:
    /**
     * Create an optimal control problem using a system, a cost, and 
     * a trajectory given by a sequence of times, states, and controls. 
     * The times ts must be given, the initial state xs[0] must be set, and
     * the optimization will be performed over the controls us
     * 
     * @param sys system
     * @param cost cost
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param p (np-size) system parameters (set to 0 if none)
     * @param dus (N) initial control standard deviations (determines how widespread the initial distro is)
     * @param es (N) zero-mean regularization/local-exploration noise variances
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    SystemCe(System<T, n, c, np> &sys, Cost<T, n, c, np, Tc> &cost, 
             vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, Vectormd *p,
             vector<Vectorcd> &dus, vector<Vectorcd> &es, 
             bool update = true);

    /**
     * Create an optimal control problem using a system, a cost, and 
     * a trajectory parametrized using a finite set of parameters of dimension ntp
     * The parametrization Tparam defines a mapping between the parameters and a 
     * discrete trajectory at times ts, states xs, and control us (and optionally system
     * parameters p)
     * Optimization occurs over the parametrization.
     *
     * 
     * @param sys system
     * @param cost cost
     * @param tp trajectory parametrization
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param p (np-size) system parameters (set to 0 if none)
     * @param dus (N) initial control standard deviations (determines how widespread the initial distro is)
     * @param es (N) zero-mean regularization/local-exploration noise variances
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    SystemCe(System<T, n, c, np> &sys, Cost<T, n, c, np, Tc> &cost, Tparam<T, n, c, np, ntp, Tc> &tp,
             vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, Vectormd *p,
             vector<Vectorcd> &dus, vector<Vectorcd> &es, 
             bool update = true);

    /**
     * Create an optimal control problem using a system, a cost, and 
     * a trajectory parametrized using a finite set of parameters of dimension ntp
     * The parametrization Tparam defines a mapping between the parameters and a 
     * discrete trajectory at times ts, states xs, and control us (and optionally system
     * parameters p)
     * Optimization occurs over the parametrization.
     *
     * 
     * @param sys system
     * @param cost cost
     * @param tp trajectory parametrization
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param p (np-size) system parameters (set to 0 if none)
     * @param dss (N) initial tparam variances (determines how widespread the initial distro is)
     * @param dess (N) zero-mean regularization/local-exploration noise variances
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    SystemCe(System<T, n, c, np> &sys, Cost<T, n, c, np, Tc> &cost, Tparam<T, n, c, np, ntp, Tc> &tp,
             vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, Vectormd *p,
             VectorXd &dss, VectorXd &dess, 
             bool update = true);

    /**
     * Create an optimal control problem using a system, a cost, and 
     * a trajectory parametrized using a finite set of parameters of dimension ntp
     * The parametrization Tparam defines a mapping between the parameters and a 
     * discrete trajectory at times ts, states xs, and control us (and optionally system
     * parameters p)
     * For this cosntructor, it is not necessary to have a mapping from discrete trajectory
     * to tparam, but only fro tparam to discrete trajectory.
     *
     * 
     * @param sys system
     * @param cost cost
     * @param tp trajectory parametrization
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param p (np-size) system parameters (set to 0 if none)
     * @param mu0 (N) initial tparam parameters mean
     * @param P0 (N) initial tparam parameters covaraiance
     * @param S (N) zero-mean regularization/local-exploration noise variances
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    SystemCe(System<T, n, c, np> &sys, Cost<T, n, c, np, Tc> &cost, Tparam<T, n, c, np, ntp, 
             Tc> &tp, Creator<Tc> *contextSampler,
             vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, Vectormd *p,
             Vectortpd &mu0, Matrixtpd &P0, Matrixtpd &S, bool update = true);

    
    
    virtual ~SystemCe();
    

    /**
     * Perform one SYSTEMCE iteration. 
     */
    void Iterate(bool updatexsfromus = true);

    /**
     * Generate a full trajectory (xs and us) from a parameter z and (optionally) return its cost
     * @param xs trajectory
     * @param us controls
     * @param evalCost whether to compute and return the trajectory cost
     */
    double Update(std::vector<T> &xs, const std::vector<Vectorcd> &us, bool evalCost = true);
        
    void us2z(Vectortpd &z, const std::vector<Vectorcd> &us) const;    

    void z2us(std::vector<Vectorcd> &us, const Vectortpd &z) const;

    System<T, n, c, np> &sys;    ///< dynamical system

    Cost<T, n, c, np, Tc> &cost;     ///< given cost function

    Tparam<T, n, c, np, ntp, Tc> *tp;     ///< trajectory parametrization

    Creator<Tc> *contextSampler;     ///< context sampler

    std::vector<double> &ts; ///< times (N+1) vector

    std::vector<T> &xs;      ///< states (N+1) vector

    std::vector<Vectorcd> &us;      ///< controls (N) vector

    Vectormd *p;      ///< system internal parameters

    std::vector<Vectorcd> &dus;      ///< control standard deviations (N) vector

    std::vector<T> xss;             ///< states (N+1) vector
    
    std::vector<Vectorcd> uss;      ///< controls (N) vector

    std::vector< std::vector<T> > xsas;    ///< all state samples
    
    std::vector< std::vector<Vectorcd> >usas;  ///< all control samples

    std::vector< bool >bs;        ///< all samples feasibility
    
    int N;        ///< number of discrete trajectory segments

    Ce<ntp> ce;       ///< GMM-based cross-entropy algorithm

    int Ns;      ///< number of cross-entropy samples

    double J;    ///< optimal cost

    bool enforceUpperBound; ///< whether to discard any samples with cost > enforceUpperBoundFactor*Jub (true by default)

    double enforceUpperBoundFactor;  ///< default is 1

    bool updateUpperBound; ///< after every iteration Jub is updated to be the current max cost (true by default)

    double Jub;  ///< allowed cost upper bound (inf by default)

    Vectortpd zmin; ///< Current best z corresponding to J

    int nofevaluations;///< Number of evaluations at any point of time
    
    bool debug;  ///< whether to display debugging info    

    // TODO: fix this: call it an external callback, etc... 
    RenderFunc *external_render;///<RenderFunction for rendering samples

  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int n, int c, int np, int ntp, typename Tc> 
    SystemCe<T, n, c, np, ntp, Tc>::SystemCe(System<T, n, c, np> &sys, 
                                             Cost<T, n, c, np, Tc> &cost, 
                                             vector<double> &ts, 
                                             vector<T> &xs, 
                                             vector<Matrix<double, c, 1> > &us,
                                             Matrix<double, np, 1> *p,
                                             vector<Matrix<double, c, 1> > &dus,
                                             vector<Matrix<double, c, 1> > &es,
                                             bool update) : 
    sys(sys), cost(cost), tp(0), contextSampler(0),
    ts(ts), xs(xs), us(us), p(p), dus(dus), xss(xs), uss(us), N(us.size()), 
    ce(N*sys.U.n, 1), Ns(1000), 
    J(numeric_limits<double>::max()), 
    enforceUpperBound(true), enforceUpperBoundFactor(1), updateUpperBound(true), Jub(numeric_limits<double>::max()),
    debug(true), external_render(0), nofevaluations(0)
    {
      // do not use an external tparam, just assume discrete controls are the params
      
      if (ntp != Dynamic) {
        assert(ntp == N*sys.U.n);
        if (ntp > 16)
          cout << "[W] gcop::SystemCe::SystemCe: fixed Eigen size is not recommended above 16!" << endl;
      }

      assert(N > 0);
      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      assert(us.size() == N);
      assert(dus.size() == N);
      assert(xss.size() == N+1);
      assert(uss.size() == N);

      if (update) {
        J = Update(xs, us);
      }
      
      us2z(ce.gmm.ns[0].mu, us);
      Vectortpd z;
      if (ntp == Dynamic)
        z.resize(sys.U.n*N);
  
      us2z(z, dus);
      ce.gmm.ns[0].P = (z.cwiseProduct()).asDiagonal();
      ce.gmm.Update();

      us2z(z, es);
      ce.S = z.asDiagonal();
    }

  template <typename T, int n, int c, int np, int ntp, typename Tc> 
    SystemCe<T, n, c, np, ntp, Tc>::SystemCe(System<T, n, c, np> &sys, 
                                             Cost<T, n, c, np, Tc> &cost, 
                                             Tparam<T, n, c, np, ntp, Tc> &tp,
                                             vector<double> &ts, 
                                             vector<T> &xs, 
                                             vector<Vectorcd> &us, 
                                             Vectormd *p,
                                             VectorXd &dss, 
                                             VectorXd &dess, 
                                             bool update) :
    sys(sys), cost(cost), tp(&tp), contextSampler(0),
    ts(ts), xs(xs), us(us), dus(us), p(p), xss(xs), uss(us), N(us.size()), 
    ce(tp.ntp, 1), Ns(1000), 
    J(numeric_limits<double>::max()), 
    enforceUpperBound(true), enforceUpperBoundFactor(1), updateUpperBound(true), Jub(numeric_limits<double>::max()),
    debug(true), external_render(0), nofevaluations(0)
    {

      assert(N > 0);
      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      assert(us.size() == N);
      assert(xss.size() == N+1);
      assert(uss.size() == N);
      assert(dss.size() == tp.ntp);
      assert(dess.size() == tp.ntp);

      if (update) {
        J = Update(xs, us);
      }
    
      tp.To(ce.gmm.ns[0].mu, ts, xs, us, p);

      ce.gmm.ns[0].P = dss.asDiagonal();
      ce.gmm.Update();

      ce.S = dess.asDiagonal();
    }

  template <typename T, int n, int c, int np, int ntp, typename Tc> 
    SystemCe<T, n, c, np, ntp, Tc>::SystemCe(System<T, n, c, np> &sys, 
                                             Cost<T, n, c, np, Tc> &cost, 
                                             Tparam<T, n, c, np, ntp, Tc> &tp,
                                             Creator<Tc> *contextSampler,
                                             vector<double> &ts, 
                                             vector<T> &xs, 
                                             vector<Vectorcd> &us, 
                                             Vectormd *p,
                                             Vectortpd &mu0,
                                             Matrixtpd &P0,
                                             Matrixtpd &S, 
                                             bool update) :
    sys(sys), cost(cost), tp(&tp), contextSampler(contextSampler),
    ts(ts), xs(xs), us(us), dus(us), p(p), xss(xs), uss(us), N(us.size()), 
    ce(tp.ntp, 1), Ns(1000), 
    J(numeric_limits<double>::max()), 
    enforceUpperBound(true), enforceUpperBoundFactor(1), updateUpperBound(true), Jub(numeric_limits<double>::max()),
    debug(true), external_render(0), nofevaluations(0)
    {

      assert(N > 0);
      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      assert(us.size() == N);
      assert(xss.size() == N+1);
      assert(uss.size() == N);
      assert(mu0.size() == tp.ntp);
      //      assert(tpCov.size() == tp.ntp);
      //      assert(tpS.size() == tp.ntp);

      //      if (update) {
      //        J = Update(xs, us);
      //      } else {
      //      }
      
      ce.gmm.ns[0].mu = mu0;
      ce.gmm.ns[0].P = P0;
      ce.S = S;

      ce.gmm.Update();
          
      //      tp.To(ce.gmm.ns[0].mu, ts, xs, us, p);
      if (contextSampler) {
        tp.SetContext( (*contextSampler)() );
      }
      tp.From(ts, xs, us, mu0, p);
    }


  template <typename T, int n, int c, int np, int ntp, typename Tc> 
    SystemCe<T, n, c, np, ntp, Tc>::SystemCe(System<T, n, c, np> &sys, 
                                         Cost<T, n, c, np, Tc> &cost,                              
                                         Tparam<T, n, c, np, ntp, Tc> &tp, 
                                         vector<double> &ts, 
                                         vector<T> &xs, 
                                         vector<Matrix<double, c, 1> > &us,
                                         Matrix<double, np, 1> *p,
                                         vector<Matrix<double, c, 1> > &dus,
                                         vector<Matrix<double, c, 1> > &es,
                                         bool update) : 
    sys(sys), cost(cost), tp(&tp), contextSampler(0),
    ts(ts), xs(xs), us(us), p(p), dus(dus), xss(xs), uss(us), N(us.size()), 
    ce(tp.ntp, 1), Ns(1000), 
    J(numeric_limits<double>::max()), 
    enforceUpperBound(true), enforceUpperBoundFactor(1), updateUpperBound(true), Jub(numeric_limits<double>::max()),
    debug(true), external_render(0), nofevaluations(0)
    {
      /*
      if (ntp != Dynamic) {
        assert(ntp == N*sys.U.n);
        if (ntp > 16)
          cout << "[W] gcop::SystemCe::SystemCe: fixed Eigen size is not recommended above 16!" << endl;
      }
      */

      assert(N > 0);
      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      assert(us.size() == N);
      assert(dus.size() == N);
      assert(xss.size() == N+1);
      assert(uss.size() == N);

      if (update) {
        J = Update(xs, us);
      }

      // initialize mean using the provided initial trajectory and controls
      tp.To(ce.gmm.ns[0].mu, ts, xs, us, p);
      Vectortpd z;
      Vectortpd dz;
      if (ntp == Dynamic) {
        z.resize(tp.ntp);
        dz.resize(tp.ntp);
      }  

      // perturbed controls:  ups = us + dus
      vector<Matrix<double, c, 1> > ups(us.size());
      transform(us.begin(), us.end(), dus.begin(), ups.begin(),plus< Matrix<double, c, 1> >());

      // initialize the covariance using the provided variance in the initial controls
      tp.To(z, ts, xs, ups, p);
      dz = z - ce.gmm.ns[0].mu;
      ce.gmm.ns[0].P = (dz.cwiseProduct(dz)).asDiagonal();
      ce.gmm.Update();

      tp.To(z, ts, xs, es, p);
      ce.S = z.asDiagonal();
    }

  
  template <typename T, int n, int c, int np, int ntp, typename Tc> 
    SystemCe<T, n, c, np, ntp, Tc>::~SystemCe()
    {
    }

  template <typename T, int n, int c, int np, int ntp, typename Tc> 
    void SystemCe<T, n, c, np, ntp, Tc>::us2z(Vectortpd &z, const std::vector<Vectorcd> &us) const
    {
      assert(us.size() > 0);
      assert(z.size() == us.size()*us[0].size());
      for (int i = 0; i < us.size(); ++i) {
        z.segment(i*sys.U.n, sys.U.n) = us[i];
      }
    }
  
  template <typename T, int n, int c, int np, int ntp, typename Tc> 
    void SystemCe<T, n, c, np, ntp, Tc>::z2us(std::vector<Vectorcd> &us, const Vectortpd &z) const
    {
      assert(us.size() > 0);
      assert(z.size() == us.size()*us[0].size());
      for (int i = 0; i < us.size(); ++i) {
        us[i] = z.segment(i*sys.U.n, sys.U.n);
      }
    }

  
  template <typename T, int n, int c, int np, int ntp, typename Tc> 
    double SystemCe<T, n, c, np, ntp, Tc>::Update(vector<T> &xs, const vector<Vectorcd> &us, bool evalCost) {    
    double J = 0;
    sys.Reset(xs[0],ts[0]);//gives a chance for physics engines to reset to specific state and time.

    for (int k = 0; k < N; ++k) {
      double h = ts[k+1] - ts[k];
      sys.Step(xs[k+1], us[k], h, p);
      if (evalCost) 
        J += cost.L(ts[k], xs[k], us[k], h, p);
      
    }
    if (evalCost)
      J += cost.L(ts[N], xs[N], us[N-1], 0, p);
    return J;
  }

  template <typename T, int n, int c, int np, int ntp, typename Tc> 
    void SystemCe<T, n, c, np, ntp, Tc>::Iterate(bool updatexsfromus) {

      if (ce.inc) 
      {
        Vectortpd z;
        if (ntp == Dynamic)
          if (tp)
            z.resize(tp->ntp);
          else 
            z.resize(us.size()*this->N);        

        ce.Sample(z);
        if (tp) {
          // set context (e.g. goal, environment, etc...) possibly using random sampling
          if (contextSampler) {
            Tc &context = (*contextSampler)();
            tp->SetContext(context);            
            cost.SetContext(context);
          }
          tp->From(ts, xss, uss, z, p);
          double cost_trajectory = 0;
          int N = uss.size();
          for(int k = 0;k < N; k++)
          { 
            cost_trajectory += cost.L(ts[k], xss[k], uss[k], ts[k+1]-ts[k], p);
          }
          cost_trajectory += cost.L(ts[N], xss[N], uss[N-1], 0, p);
          ce.AddSample(z, cost_trajectory);
        }
        else
        {
          z2us(uss, z);
          ce.AddSample(z, Update(xss, uss));
        }

        xsas.push_back(xss);
        usas.push_back(uss);
      } else   // non-incremental
      {
        // Jub is used to discard any samples with cost below Jub
        ce.Reset();
        xsas.clear();
        usas.clear();
        bs.clear();

        Vectortpd z;
        if (ntp == Dynamic)
          if (tp)
            z.resize(tp->ntp);
          else 
            z.resize(us.size()*this->N);
        
        for (int j = 0; j < Ns; ++j) {
          
          // sample context
          if (contextSampler) {
            Tc &context = (*contextSampler)();
            tp->SetContext(context);            
            cost.SetContext(context);
          }
          
          double J = 0;
          bool good = true;

          do {
            ce.Sample(z);       
            if (tp) {
              good = tp->From(ts, xss, uss, z, p);
              int N = uss.size();
              J = 0;
              for(int k = 0;k < N; k++) { 
                J += cost.L(ts[k], xss[k], uss[k], ts[k+1]-ts[k], p);
              }
              J += cost.L(ts[N], xss[N], uss[N-1], 0, p);
            } else {
              z2us(uss, z);
              J = Update(xss, uss);
            }
          } while(enforceUpperBound && J > enforceUpperBoundFactor*Jub);
          
          ce.AddSample(z, J);
          
          // add to list of all samples
          xsas.push_back(xss);
          usas.push_back(uss);
          bs.push_back(good);

          //#DEBUG Number of function evaluations:
          //cout<<++nofevaluations<<endl;
          ++nofevaluations;

          //Render trajectory samples if external rendering function is provided:
          if(external_render)
          {
            external_render(j,xss);//ID for the sample trajectory
          }
        }
        
        if (updateUpperBound) {
          //cout << "Jub=" << Jub << " ce.Jmax=" << ce.Jmax << endl;          
          Jub = min(Jub, ce.Jmax);
          //cout << "Jub=" << Jub << endl;
        }
      }

    // estimate distribution
    ce.Select();    

    if (!ce.Fit()) {
      cout << "[W] TrajectoryPrmSample::Sample: ce.Fit failed!" << endl;
    }

    if (ce.Jmin < J) {
      if (tp)
        tp->From(ts, xs, us, ce.zmin, p);
      else
        z2us(us, ce.zmin);

      if(updatexsfromus)
        Update(xs, us, false);
      J = ce.Jmin;
      zmin = ce.zmin;//Store the parameters also
    }

    // construct trajectory using the first sample (this is the one with lowest cost)
    //    Traj(xs, ce.zps[0].first, x0);
    // cout << "Solution: xs[K]=" << xs[K].transpose() << " c=" << ce.cs[0] << endl;
  }
}

#endif
