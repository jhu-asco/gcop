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
 
  template <typename T, int n = Dynamic, int c = Dynamic> class Docp {
    
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
    Docp(System<T, Vectorcd, n, c> &sys, Cost<T, Vectorcd, n, c> &cost, 
         vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, bool update = true);
    
    virtual ~Docp();


    /**
     * Perform one DOCP iteration. Internally calls:
     *
     * Backward -> Forward -> Update. The controls us and trajectory xs
     * are updated. 
     */
    void Iterate();

    /**
     * Update the trajectory and (optionally) its linearization
     * @param der whether to update derivatives (A and B matrices)
     */
    void Update(bool der = true);    
    
    System<T, Vectorcd, n, c> &sys;    ///< dynamical system 

    Cost<T, Vectorcd, n, c> &cost;     ///< given cost function

    std::vector<double> &ts; ///< times (N+1) vector

    std::vector<T> &xs;      ///< states (N+1) vector

    std::vector<Vectorcd> &us;      ///< controls (N) vector
    
    int N;        ///< number of discrete trajectory segments
    
    double mu;    ///< current regularization factor mu
    double mu0;   ///< minimum regularization factor mu
    double dmu0;  ///< regularization factor modification step-size
    double mumax; ///< maximum regularization factor mu
    double a;     ///< step-size
    
    std::vector<Vectorcd> dus;  ///< computed control change
    
    std::vector<Matrixnd> As;
    std::vector<Matrixncd> Bs;
    
    std::vector<Vectorcd> kus;
    std::vector<Matrixcnd> Kuxs;
    
    Vectornd Lx;
    Matrixnd Lxx;
    Vectorcd Lu;
    Matrixcd Luu;
    
    Matrixnd P;
    Vectornd v;
    
    double V;
    Vector2d dV;

    double s1;   ///< Armijo/Bertsekas step-size control factor s1
    double s2;   ///< Armijo/Bertsekas step-size control factor s2
    double b1;   ///< Armijo/Bertsekas step-size control factor b1
    double b2;   ///< Armijo/Bertsekas step-size control factor b2

    bool debug;  ///< whether to display debugging info
    
    double eps;                     ///< epsilon used for finite differences

    bool pd(const Matrixnd &P) {
      LLT<Matrixnd> llt;     
      llt.compute(P);
      return llt.info() == Eigen::Success;
    }

    bool pdX(const MatrixXd &P) {
      LLT<MatrixXd> llt;     
      llt.compute(P);
      return llt.info() == Eigen::Success;
    }  


  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int n, int c> 
    Docp<T, n, c>::Docp(System<T, Matrix<double, c, 1>, n, c> &sys, 
                        Cost<T, Matrix<double, c, 1>, n, c> &cost, 
                        vector<double> &ts, 
                        vector<T> &xs, 
                        vector<Matrix<double, c, 1> > &us,
                        bool update) : 
    sys(sys), cost(cost), ts(ts), xs(xs), us(us), N(us.size()), 
    mu(1e-3), mu0(1e-3), dmu0(2), mumax(1e6), a(1), 
    dus(N),
    As(N), Bs(N), kus(N), Kuxs(N), 
    s1(0.1), s2(0.5), b1(0.25), b2(2),
    debug(true), eps(1e-6)
    {
      assert(N > 0);
      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      assert(us.size() == N);

      if (n == Dynamic || c == Dynamic) {
        for (int i = 0; i < N; ++i) {
          dus[i].resize(sys.c);
          As[i].resize(sys.n, sys.n);
          Bs[i].resize(sys.n, sys.c);
          kus[i].resize(sys.c);
          Kuxs[i].resize(sys.c, sys.n);
        }

        Lx.resize(sys.n);
        Lxx.resize(sys.n, sys.n);
        Lu.resize(sys.c);
        Luu.resize(sys.c, sys.c);

        P.resize(sys.n, sys.n);
        v.resize(sys.n);       
      }

      if (update) {
        Update();
      }
    }
  
  template <typename T, int n, int c> 
    Docp<T, n, c>::~Docp()
    {
    }
  
  template <typename T, int n, int c> 
    void Docp<T, n, c>::Update(bool der) {

    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    
    Vectornd dx;
    Vectorcd du;
    Vectornd dfp;
    Vectornd dfm;
    T xav(xs[0]);
    T xbv(xs[0]);
    
    //    cout << "SIZE=" << ((MbsState*)&xav)->r.size() << endl;
    
    if (n == Dynamic) {
      dx.resize(sys.n);
      dfp.resize(sys.n);
      dfm.resize(sys.n);
    }
    
    if (c == Dynamic)
      du.resize(sys.c);    

    for (int k = 0; k < N; ++k) {
      double h = ts[k+1] - ts[k];
      if (der) {
        static const double q = 1093121342312;  // a random number
        As[k](0,0) = q;
        Bs[k](0,0) = q;
        
        const T &xa = xs[k];
        T &xb = xs[k+1];
        const Vectorcd &u = us[k];

        sys.Step(xb, ts[k], xa, u, h, 0, &As[k], &Bs[k], 0);

        //        cout << "B=" << endl << Bs[k] << endl;
        

        //        if (fabs(As[k](0,0) - q) < 1e-10) {
        //          autodiff.DF(xb, ts[k], xa, u, h, &As[k], &Bs[k]);
        //        }

        
        // if no jacobians were provided use finite differences
        if (fabs(As[k](0,0) - q) < 1e-10) {
                    
          assert(sys.n > 0 && sys.c > 0);

          for (int i = 0; i < sys.n; ++i) {
            dx.setZero();
            dx[i] = eps;
            
            // xav = xa + dx
            cost.X.Retract(xav, xa, dx);
            
            // reconstruct state using previous time-step
            sys.Rec(xav, h);
            
            // xbv = F(xav)
            sys.Step(xbv, ts[k], xav, u, h);
            
            // df = xbv - xb
            cost.X.Lift(dfp, xb, xbv);

            // xav = xa - dx
            cost.X.Retract(xav, xa, -dx);
            
            // reconstruct state using previous time-step
            sys.Rec(xav, h);
            
            // xbv = F(xav)
            sys.Step(xbv, ts[k], xav, u, h);
            
            // dfm = xbv - xb
            cost.X.Lift(dfm, xb, xbv);

            As[k].col(i) = (dfp - dfm)/(2*eps);
          }
        }
        
        if (fabs(Bs[k](0,0) - q) < 1e-10) {
          
          for (int i = 0; i < sys.c; ++i) {
            du.setZero();
            du[i] = eps;
            
            // xbv = F(xa, u + du)
            sys.Step(xbv, ts[k], xa, u + du, h);
            
            // df = xbv - xb
            cost.X.Lift(dfp, xb, xbv);

            // xbv = F(xa, u - du)
            sys.Step(xbv, ts[k], xa, u - du, h);
            
            // df = xbv - xb
            cost.X.Lift(dfm, xb, xbv);

            Bs[k].col(i) = (dfp - dfm)/eps;
          }
        }
        
      } else {
        sys.Step(xs[k+1], ts[k], xs[k], us[k], h);
      }
    }
  }  

  template <typename T, int n, int c> 
    void Docp<T, n, c>::Iterate() {
  }
}

#endif
