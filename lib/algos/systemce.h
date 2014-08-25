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

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
  template <typename T, int n = Dynamic, int c = Dynamic> class SystemCe {
    
    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, n, n> Matrixnd;
    typedef Matrix<double, n, c> Matrixncd;
    typedef Matrix<double, c, n> Matrixcnd;
    typedef Matrix<double, c, c> Matrixcd;  
    
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
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    SystemCe(System<T, n, c> &sys, Cost<T, n, c> &cost, 
             vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
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
        
    void us2z(VectorXd &z, const std::vector<Vectorcd> &us) const;    

    void z2us(std::vector<Vectorcd> &us, const VectorXd &z) const;

    System<T, n, c> &sys;    ///< dynamical system

    Cost<T, n, c> &cost;     ///< given cost function

    std::vector<double> &ts; ///< times (N+1) vector

    std::vector<T> &xs;      ///< states (N+1) vector

    std::vector<Vectorcd> &us;      ///< controls (N) vector

    std::vector<Vectorcd> &dus;      ///< control variation (N) vector

    std::vector<T> xss;             ///< states (N+1) vector
    
    std::vector<Vectorcd> uss;      ///< controls (N) vector
    
    int N;        ///< number of discrete trajectory segments

    Ce ce;       ///< GMM-based cross-entropy algorithm

    int Ns;      ///< number of cross-entropy samples

    double J;    ///< optimal cost
    
    bool debug;  ///< whether to display debugging info    

  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int n, int c> 
    SystemCe<T, n, c>::SystemCe(System<T, n, c> &sys, 
                                Cost<T, n, c> &cost, 
                                vector<double> &ts, 
                                vector<T> &xs, 
                                vector<Matrix<double, c, 1> > &us,
                                vector<Matrix<double, c, 1> > &dus,
                                vector<Matrix<double, c, 1> > &es,
                                bool update) : 
    sys(sys), cost(cost), ts(ts), xs(xs), us(us), dus(dus), xss(xs), uss(us), N(us.size()), 
    ce(N*sys.U.n, 1), Ns(1000), debug(true)
    {
      assert(N > 0);
      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      assert(us.size() == N);
      assert(dus.size() == N);
      assert(xss.size() == N+1);
      assert(uss.size() == N);
      
      if (n == Dynamic || c == Dynamic) {
        /*
        for (int i = 0; i < N; ++i) {
          dus[i].resize(sys.U.n);
          As[i].resize(sys.n, sys.n);
          Bs[i].resize(sys.n, sys.U.n);
          kus[i].resize(sys.U.n);
          Kuxs[i].resize(sys.U.n, sys.n);
        }
        */
      }

      if (update) {
        J = Update(xs, us);
      } else {
        J = numeric_limits<double>::max();
      }
      
      us2z(ce.gmm.ns[0].mu, us);
      VectorXd zdu(sys.U.n*N);
      us2z(zdu, dus);
      ce.gmm.ns[0].P = zdu.asDiagonal();
      ce.gmm.Update();
      VectorXd ze(sys.U.n*N);
      us2z(ze, es);
      ce.S = ze.asDiagonal();
    }
  
  template <typename T, int n, int c> 
    SystemCe<T, n, c>::~SystemCe()
    {
    }

  template <typename T, int n, int c> 
    void SystemCe<T, n, c>::us2z(VectorXd &z, const std::vector<Vectorcd> &us) const
    {
      assert(us.size() > 0);
      assert(z.size() == us.size()*us[0].size());
      for (int i = 0; i < us.size(); ++i) {
        z.segment(i*sys.U.n, sys.U.n) = us[i];
      }
    }
  
  template <typename T, int n, int c> 
    void SystemCe<T, n, c>::z2us(std::vector<Vectorcd> &us, const VectorXd &z) const
    {
      assert(us.size() > 0);
      assert(z.size() == us.size()*us[0].size());
      for (int i = 0; i < us.size(); ++i) {
        us[i] = z.segment(i*sys.U.n, sys.U.n);
      }
    }

  
  template <typename T, int n, int c> 
    double SystemCe<T, n, c>::Update(vector<T> &xs, const vector<Vectorcd> &us, bool evalCost) {    
    double J = 0;
    // @mk:TODO fix that
    // sys.reset();//gives a chance for physics engines to reset themselves. Added Gowtham 8/2/14

    for (int k = 0; k < N; ++k) {
      double h = ts[k+1] - ts[k];
      sys.Step(xs[k+1], ts[k], xs[k], us[k], h);
      if (evalCost) 
        J += cost.L(ts[k], xs[k], us[k], h);
      
    }
    if (evalCost)
      J += cost.L(ts[N], xs[N], us[N-1], 0);
    return J;
  }

  template <typename T, int n, int c> 
    void SystemCe<T, n, c>::Iterate() {

    ce.Reset();
    
    VectorXd z(sys.U.n*N);            // parameter vector

    for (int j = 0; j < Ns; ++j) {
      ce.Sample(z);
      z2us(uss, z);
      ce.AddSample(z, Update(xss, uss));
    }

    // estimate distribution
    ce.Select();    
    if (!ce.Fit()) {
      cout << "[W] TrajectoryPrmSample::Sample: ce.Fit failed!" << endl;
    }
    
    if (ce.cs[0] < J) {
      z2us(us, ce.zps[0].first);
      Update(xs, us, false);
      J = ce.cs[0];
    }

    // construct trajectory using the first sample (this is the one with lowest cost)
    //    Traj(xs, ce.zps[0].first, x0);
    // cout << "Solution: xs[K]=" << xs[K].transpose() << " c=" << ce.cs[0] << endl;
  }
}

#endif
