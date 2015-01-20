#ifndef GCOP_DOCP_H
#define GCOP_DOCP_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include "system.h"
#include "cost.h"
#include <cmath>
#include "rn.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic> class Docp {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  
    
  public:
    /**
     * Create a discrete optimal control problem using a system, a cost, and 
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
    Docp(System<T, nx, nu, np> &sys, Cost<T, nx, nu, np> &cost, 
         vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
         Vectormd *p = 0, bool update = true);
    
    virtual ~Docp();

    /**
     * Perform one DOCP iteration. 
     */
    virtual void Iterate();

    /**
     * Update the trajectory and (optionally) its linearization
     * @param der whether to update derivatives (A and B matrices)
     */
    void Update(bool der = true);
    
    System<T, nx, nu, np> &sys;    ///< dynamical system 

    Cost<T, nx, nu, np> &cost;     ///< given cost function

    std::vector<double> &ts; ///< times (N+1) vector

    std::vector<T> &xs;      ///< states (N+1) vector

    std::vector<Vectorcd> &us;      ///< controls (N) vector
    
    Vectormd *p;               ///< parameter vector
    
    std::vector<Matrixnd> As;    ///< state jacobians along the path (computed internally)
    std::vector<Matrixncd> Bs;   ///< control jacobians along the path (computed internally)
    
    bool debug;  ///< whether to display debugging info
    
    double eps;  ///< epsilon used for finite differences

    double J;    ///< Optimal cost at any point of time

    int nofevaluations;///< Number of Function evaluations at any point of time

  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int nx, int nu, int np> 
    Docp<T, nx, nu, np>::Docp(System<T, nx, nu, np> &sys, 
                              Cost<T, nx, nu, np> &cost, 
                              vector<double> &ts, 
                              vector<T> &xs, 
                              vector<Matrix<double, nu, 1> > &us,
                              Matrix<double, np, 1> *p,
                              bool update) : 
    sys(sys), cost(cost), ts(ts), xs(xs), us(us), p(p),
    As(us.size()), Bs(us.size()), J(std::numeric_limits<double>::max()), 
    debug(true), eps(1e-3), nofevaluations(0)
    {
      int N = us.size();
      assert(N > 0);
      assert(ts.size() == N+1);
      assert(xs.size() == N+1);

      if (nx == Dynamic || nu == Dynamic) {
        for (int i = 0; i < N; ++i) {
          As[i].resize(sys.X.n, sys.X.n);
          Bs[i].resize(sys.X.n, sys.U.n);
        }
      }

      if (update) {
        Update();
      }
    }
  
  template <typename T, int nx, int nu, int np> 
    Docp<T, nx, nu, np>::~Docp()
    {
    }
  
  template <typename T, int nx, int nu, int np> 
    void Docp<T, nx, nu, np>::Update(bool der) {

    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    
    Vectornd dx;
    Vectorcd du;
    Vectornd dfp;
    Vectornd dfm;
    T xav(xs[0]);
    T xbv(xs[0]);
    
    //    cout << "SIZE=" << ((MbsState*)&xav)->r.size() << endl;
    
    if (nx == Dynamic) {
      dx.resize(sys.X.n);
      dfp.resize(sys.X.n);
      dfm.resize(sys.X.n);
    }
    
    if (nu == Dynamic)
      du.resize(sys.U.n);    

    int N = us.size();

    sys.reset(xs[0],ts[0]);//Reset the internal state

    for (int k = 0; k < N; ++k) {
      double h = ts[k+1] - ts[k];
      if (der) {
        static const double q = 1093121342312;  // a random number
        As[k](0,0) = q;
        Bs[k](0,0) = q;
        
        const T &xa = xs[k];
        T &xb = xs[k+1];
        const Vectorcd &u = us[k];

        //sys.Step(xb, ts[k], xa, u, h, p, &As[k], &Bs[k], 0);
        sys.Step1(xs[k+1], us[k], h, p, &As[k], &Bs[k], 0);

        //        cout << "B=" << endl << Bs[k] << endl;
        

        //        if (fabs(As[k](0,0) - q) < 1e-10) {
        //          autodiff.DF(xb, ts[k], xa, u, h, &As[k], &Bs[k]);
        //        }

        
        // if no jacobians were provided use finite differences
        if (fabs(As[k](0,0) - q) < 1e-10) {
                    
          assert(sys.X.n > 0 && sys.U.n > 0);

          for (int i = 0; i < sys.X.n; ++i) {
            dx.setZero();
            dx[i] = eps;
            
            // xav = xa + dx
            sys.X.Retract(xav, xa, dx);
            
            // reconstruct state using previous time-step
            sys.Rec(xav, h);
            
            // xbv = F(xav)
            sys.Step(xbv, ts[k], xav, u, h, p);
            
            // df = xbv - xb
            sys.X.Lift(dfp, xb, xbv);

            // xav = xa - dx
            sys.X.Retract(xav, xa, -dx);
            
            // reconstruct state using previous time-step
            sys.Rec(xav, h);
            
            // xbv = F(xav)
            sys.Step(xbv, ts[k], xav, u, h, p);
            
            // dfm = xbv - xb
            sys.X.Lift(dfm, xb, xbv);

            As[k].col(i) = (dfp - dfm)/(2*eps);
          }
          sys.reset(xs[k+1],ts[k+1]);//Reset the internal state
          //cout<<"As["<<k<<"]: "<<endl<<As[k]<<endl;
        }
        
        if (fabs(Bs[k](0,0) - q) < 1e-10) {
          
          for (int i = 0; i < sys.U.n; ++i) {
            du.setZero();
            du[i] = eps;
            
            // xbv = F(xa, u + du)
            sys.Step(xbv, ts[k], xa, u + du, h, p);
            
            // df = xbv - xb
            sys.X.Lift(dfp, xb, xbv);

            // xbv = F(xa, u - du)
            sys.Step(xbv, ts[k], xa, u - du, h, p);
            
            // df = xbv - xb
            sys.X.Lift(dfm, xb, xbv);

            Bs[k].col(i) = (dfp - dfm)/eps;
          }
          sys.reset(xs[k+1],ts[k+1]);//Reset the internal state
         // cout<<"Bs["<<k<<"]: "<<endl<<Bs[k]<<endl;
        }        
      } else {
        sys.Step1(xs[k+1], us[k], h, p);
      }
    }
  } 

  template <typename T, int nx, int nu, int np> 
    void Docp<T, nx, nu, np>::Iterate() {
    cout << "[W] Docp::Iterate: subclasses should implement this!" << endl;
  }
}

#endif
