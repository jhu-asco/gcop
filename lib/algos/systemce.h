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

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
  template <typename T, int n = Dynamic, int c = Dynamic, int np = Dynamic, int ntp = Dynamic> class SystemCe {
    
    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, n, n> Matrixnd;
    typedef Matrix<double, n, c> Matrixncd;
    typedef Matrix<double, c, n> Matrixcnd;
    typedef Matrix<double, c, c> Matrixcd;  

    typedef Matrix<double, ntp, 1> Vectortpd;  
    
    //External render Function for trajectories:
    typedef void(RenderFunc)(int, vector<T>&);
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
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param p (np-size) system parameters (set to 0 if none)
     * @param dus (N) initial control variations (determines how widespread the initial distro is)
     * @param es (N) zero-mean regularization/local-exploration noise variances
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    SystemCe(System<T, n, c, np> &sys, Cost<T, n, c, np> &cost, 
             vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, Vectormd *p,
             vector<Vectorcd> &dus, vector<Vectorcd> &es, 
             bool update = true);

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
     * @param tp trajectory parametrization
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param p (np-size) system parameters (set to 0 if none)
     * @param dus (N) initial control variations (determines how widespread the initial distro is)
     * @param es (N) zero-mean regularization/local-exploration noise variances
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    SystemCe(System<T, n, c, np> &sys, Cost<T, n, c, np> &cost, Tparam<T, n, c, np, ntp> &tp,
             vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, Vectormd *p,
             vector<Vectorcd> &dus, vector<Vectorcd> &es, 
             bool update = true);

    
    virtual ~SystemCe();
    

    /**
     * Perform one SYSTEMCE iteration. 
     */
    void Iterate();

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

    Cost<T, n, c, np> &cost;     ///< given cost function

    Tparam<T, n, c, np, ntp> *tp;     ///< trajectory parametrization

    std::vector<double> &ts; ///< times (N+1) vector

    std::vector<T> &xs;      ///< states (N+1) vector

    std::vector<Vectorcd> &us;      ///< controls (N) vector

    Vectormd *p;      ///< controls (N) vector

    std::vector<Vectorcd> &dus;      ///< control variation (N) vector

    std::vector<T> xss;             ///< states (N+1) vector
    
    std::vector<Vectorcd> uss;      ///< controls (N) vector
    
    int N;        ///< number of discrete trajectory segments

    Ce<ntp> ce;       ///< GMM-based cross-entropy algorithm

    int Ns;      ///< number of cross-entropy samples

    double J;    ///< optimal cost

    Vectortpd zmin; ///< Current best z corresponding to J

    int nofevaluations;///< Number of evaluations at any point of time
    
    bool debug;  ///< whether to display debugging info    

    RenderFunc *external_render;///<RenderFunction for rendering samples

  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int n, int c, int np, int ntp> 
    SystemCe<T, n, c, np, ntp>::SystemCe(System<T, n, c, np> &sys, 
                                     Cost<T, n, c, np> &cost, 
                                     vector<double> &ts, 
                                     vector<T> &xs, 
                                     vector<Matrix<double, c, 1> > &us,
                                     Matrix<double, np, 1> *p,
                                     vector<Matrix<double, c, 1> > &dus,
                                     vector<Matrix<double, c, 1> > &es,
                                     bool update) : 
    sys(sys), cost(cost), tp(0),
    ts(ts), xs(xs), us(us), p(p), dus(dus), xss(xs), uss(us), N(us.size()), 
    ce(N*sys.U.n, 1), Ns(1000), debug(true), external_render(0), nofevaluations(0)
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
      } else {
        J = numeric_limits<double>::max();
      }
      
      us2z(ce.gmm.ns[0].mu, us);
      Vectortpd z;
      if (ntp == Dynamic)
        z.resize(sys.U.n*N);
  
      us2z(z, dus);
      ce.gmm.ns[0].P = z.asDiagonal();
      ce.gmm.Update();

      us2z(z, es);
      ce.S = z.asDiagonal();
    }

  template <typename T, int n, int c, int np, int ntp> 
    SystemCe<T, n, c, np, ntp>::SystemCe(System<T, n, c, np> &sys, 
                                         Cost<T, n, c, np> &cost,                              
                                         Tparam<T, n, c, np, ntp> &tp, 
                                         vector<double> &ts, 
                                         vector<T> &xs, 
                                         vector<Matrix<double, c, 1> > &us,
                                         Matrix<double, np, 1> *p,
                                         vector<Matrix<double, c, 1> > &dus,
                                         vector<Matrix<double, c, 1> > &es,
                                         bool update) : 
    sys(sys), cost(cost), tp(&tp), 
    ts(ts), xs(xs), us(us), p(p), dus(dus), xss(xs), uss(us), N(us.size()), 
    ce(tp.ntp, 1), Ns(1000), debug(true), external_render(0), nofevaluations(0)
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
      } else {
        J = numeric_limits<double>::max();
      }
    
      tp.To(ce.gmm.ns[0].mu, ts, xs, us, p);
            //      us2z(ce.gmm.ns[0].mu, us);
      Vectortpd z;
      if (ntp == Dynamic) {
        z.resize(tp.ntp);
      }  


      tp.To(z, ts, xs, dus, p);
      // us2z(z, dus);
      ce.gmm.ns[0].P = z.asDiagonal();
      ce.gmm.Update();

      //      us2z(z, es);
      tp.To(z, ts, xs, es, p);
      ce.S = z.asDiagonal();
    }

  
  template <typename T, int n, int c, int np, int ntp> 
    SystemCe<T, n, c, np, ntp>::~SystemCe()
    {
    }

  template <typename T, int n, int c, int np, int ntp> 
    void SystemCe<T, n, c, np, ntp>::us2z(Vectortpd &z, const std::vector<Vectorcd> &us) const
    {
      assert(us.size() > 0);
      assert(z.size() == us.size()*us[0].size());
      for (int i = 0; i < us.size(); ++i) {
        z.segment(i*sys.U.n, sys.U.n) = us[i];
      }
    }
  
  template <typename T, int n, int c, int np, int ntp> 
    void SystemCe<T, n, c, np, ntp>::z2us(std::vector<Vectorcd> &us, const Vectortpd &z) const
    {
      assert(us.size() > 0);
      assert(z.size() == us.size()*us[0].size());
      for (int i = 0; i < us.size(); ++i) {
        us[i] = z.segment(i*sys.U.n, sys.U.n);
      }
    }

  
  template <typename T, int n, int c, int np, int ntp> 
    double SystemCe<T, n, c, np, ntp>::Update(vector<T> &xs, const vector<Vectorcd> &us, bool evalCost) {    
    double J = 0;
    sys.reset(xs[0],ts[0]);//gives a chance for physics engines to reset to specific state and time.

    for (int k = 0; k < N; ++k) {
      double h = ts[k+1] - ts[k];
      sys.Step_internalinput(xs[k+1], us[k], h, p);
      if (evalCost) 
        J += cost.L(ts[k], xs[k], us[k], h, p);
      
    }
    if (evalCost)
      J += cost.L(ts[N], xs[N], us[N-1], 0, p);
    return J;
  }

  template <typename T, int n, int c, int np, int ntp> 
    void SystemCe<T, n, c, np, ntp>::Iterate() {

      if (ce.inc) 
      {
        Vectortpd z;
        if (ntp == Dynamic)
          if (tp)
            z.resize(tp->ntp);
          else 
            z.resize(us.size()*this->N);        

        ce.Sample(z);
        if (tp)
        {
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
      } 
      else 
      {
        ce.Reset();

        Vectortpd z;
        if (ntp == Dynamic)
          if (tp)
            z.resize(tp->ntp);
          else 
            z.resize(us.size()*this->N);        

        for (int j = 0; j < Ns; ++j) {
          ce.Sample(z);
          if (tp)
          {
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

          //#DEBUG Number of function evaluations:
          //cout<<++nofevaluations<<endl;
          ++nofevaluations;

          //Render trajectory samples if external rendering function is provided:
          if(external_render)
          {
            external_render(j,xss);//ID for the sample trajectory
          }
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
