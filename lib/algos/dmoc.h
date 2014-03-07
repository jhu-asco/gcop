#ifndef GCOP_DMOC_H
#define GCOP_DMOC_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include "system.h"
#include "cost.h"
#include "mbsstate.h"
#include <cmath>
#include "rn.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
  template <typename T, int n = Dynamic, int c = Dynamic> class Dmoc {
    
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
    Dmoc(System<T, Vectorcd, n, c> &sys, Cost<T, Vectorcd, n, c> &cost, 
         vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, bool update = true);
    
    virtual ~Dmoc();


    /**
     * Perform one DMOC iteration. Internally calls:
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
    
    /**
     *  Forward pass
     */
    void Forward();

    /**
     *  Backward pass
     */
    void Backward();
    
    System<T, Vectorcd, n, c> &sys;    ///< dynamical system 

    Cost<T, Vectorcd, n, c> &cost;     ///< given cost function

    std::vector<double> &ts; ///< times (N+1) vector

    std::vector<T> &xs;      ///< states (N+1) vector

    std::vector<Vectorcd> &us;      ///< controls (N) vector
    
    int N;        ///< number of discrete trajectory segments
    
    double mu;    ///< current regularization factor mu
    double mu0;   ///< minimum regularization factor mu
    double dmu0;  ///< regularization factor modification step-size
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
    
    char type;   ///< type of algorithm (choices are PURE, DDP, LQS), LQS is default. In the current implementation second-order expansion of the dynamics is ignored, so DDP and LQS are identical. Both LQS and DDP are implemented using the true nonlinear dynamics integration in the Forward step, while PURE uses the linearized version in order to match exactly the true Newton step. 
    
    static const char PURE = 0;     ///< PURE version of algorithm (i.e. stage-wise Newton)
    static const char DDP = 1;      ///< DDP version of algorithm
    static const char LQS = 2;      ///< Linear-Quadratic Subproblem version of algorithm due to Dreyfus / Dunn&Bertsekas / Pantoja

    /*
    class Fx : public Function<T, n, n> {
    public:
    Fx(System &sys) :
      Function<T>(sys.cspace), sys(sys) {
      }
      
      void F(Vectornd &f, const T &x) {

        
        sys.X.L
        this->xa = xa;
        sys.FK(this->xa);
        sys.ID(f, t, *xb, this->xa, *u, h);
      }
      
      Mbs &sys;
      
      double t;
      const MbsState *xb;
      const Vectorcd *u;
      double h;
    };

    class Fu : public Function<VectorXd> {
    public:
    Fu(Mbs &sys) :
      Function<VectorXd>(sys.U), sys(sys) {
      }
      
      void F(VectorXd &f, const VectorXd &u) {
        sys.ID(f, t, *xb, *xa, u, h);
      }
      
      Mbs &sys;
      
      double t;
      const MbsState *xa;
      const MbsState *xb;
      double h;
    };

    */

  };
  
  
  using namespace std;
  using namespace Eigen;
  
  template <typename T, int n, int c> 
    Dmoc<T, n, c>::Dmoc(System<T, Matrix<double, c, 1>, n, c> &sys, 
                        Cost<T, Matrix<double, c, 1>, n, c> &cost, 
                        vector<double> &ts, 
                        vector<T> &xs, 
                        vector<Matrix<double, c, 1> > &us,
                        bool update) : 
    sys(sys), cost(cost), ts(ts), xs(xs), us(us), N(us.size()), 
    mu(1e-3), mu0(1e-3), dmu0(2), a(1), 
    dus(N),
    As(N), Bs(N), kus(N), Kuxs(N), 
    s1(0.1), s2(0.5), b1(0.25), b2(2),
    debug(true), type(LQS)
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
    Dmoc<T, n, c>::~Dmoc()
    {
    }
  
  template <typename T, int n, int c> 
    void Dmoc<T, n, c>::Update(bool der) {

    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    
    Vectornd dx;
    Vectorcd du;
    Vectornd df;
    T xav(xs[0]);
    T xbv(xs[0]);
    
    //    cout << "SIZE=" << ((MbsState*)&xav)->r.size() << endl;
    
    if (n == Dynamic) {
      dx.resize(sys.n);
      df.resize(sys.n);
    }
    
    if (c == Dynamic)
      du.resize(sys.c);    

    double eps = 1e-6;    

    for (int k = 0; k < N; ++k) {
      double h = ts[k+1] - ts[k];
      if (der) {
        static const double q = 1093121342312;  // a random number
        As[k](0,0) = q;
        Bs[k](0,0) = q;
        
        const T &xa = xs[k];
        T &xb = xs[k+1];
        const Vectorcd &u = us[k];

        sys.Step(xb, ts[k], xa, u, h, &As[k], &Bs[k]);


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
            cost.X.Lift(df, xb, xbv);
            As[k].col(i) = df/eps;
          }
        }
        
        if (fabs(Bs[k](0,0) - q) < 1e-10) {
          
          for (int i = 0; i < sys.c; ++i) {
            du.setZero();
            du[i] = eps;
            
            // xbv = F(xa, u + du)
            sys.Step(xbv, ts[k], xa, u + du, h);
            
            // df = xbv - xb
            cost.X.Lift(df, xb, xbv);
            Bs[k].col(i) = df/eps;
          }
        }
        
      } else {
        sys.Step(xs[k+1], ts[k], xs[k], us[k], h);
      }
    }
  }
  
  template <typename T, int n, int c> 
    void Dmoc<T, n, c>::Backward() {

    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, n, n> Matrixnd;
    typedef Matrix<double, n, c> Matrixncd;
    typedef Matrix<double, c, n> Matrixcnd;
    typedef Matrix<double, c, c> Matrixcd;  
    
    double t = ts.back();
    const T &x = xs.back();
    const Vectorcd &u = us.back();
    double L = cost.L(t, x, u, 0, &Lx, &Lxx);
    
    V = L;
    dV.setZero();
    
    v = Lx;
    P = Lxx;
    
    Vectornd Qx;
    Vectorcd Qu;    
    Matrixnd Qxx;
    Matrixcd Quu;
    Matrixcd Quum;
    Matrixcnd Qux;   
    
    Matrixnd At;
    Matrixcnd Bt;
    
    Matrixcd Ic = MatrixXd::Identity(sys.c, sys.c);
    
    for (int k = N-1; k >=0; --k) {
      
      t = ts[k];
      double h = ts[k+1] - t;

      const T &x = xs[k];
      const Vectorcd &u = us[k];
      double L = cost.L(t, x, u, h, &Lx, &Lxx, &Lu, &Luu);

      if (std::isnan(L)) {
        cout << "NAN" << " k=" << k << " Lx=" << Lx << " Lxx=" << Lxx << endl;
        getchar();
      }
      
      V += L;


      
      const Matrixnd &A = As[k];
      const Matrixncd &B = Bs[k];
      Vectorcd &ku = kus[k];
      Matrixcnd &Kux = Kuxs[k];
      
      At = A.transpose();    
      Bt = B.transpose();
      
      Qx = Lx + At*v;
      Qu = Lu + Bt*v;
      Qxx = Lxx + At*P*A;
      Quu = Luu + Bt*P*B;
      Qux = Bt*P*A;
      
      double mu = this->mu;
      double dmu = 1;
      
      LLT<Matrixcd> llt;
      
      while (1) {
        Quum = Quu + mu*Ic;
        llt.compute(Quum);
        
        // if OK, then reduce mu and break
        if (llt.info() == Eigen::Success) {
          // this is the quadratic rule proposed by Tassa and Todorov
          dmu = min(1/this->dmu0, dmu/this->dmu0);
          if (mu*dmu > this->mu0)
            mu = mu*dmu;
          else
            mu = this->mu0;

          break;
        }
        
        // if negative then increase mu
        dmu = max(this->dmu0, dmu*this->dmu0);
        mu = max(this->mu0, mu*dmu);   
        
        if (debug) {
          cout << "[I] Dmoc::Backward: increased mu=" << mu << " at k=" << k << endl;
          cout << "[I] Dmoc::Backward: Quu=" << Quu << endl;
        }
      }
      
      ku = -llt.solve(Qu);
      Kux = -llt.solve(Qux);
      
      assert(!std::isnan(ku[0]));
      assert(!std::isnan(Kux(0,0)));

      v = Qx + Kux.transpose()*Qu;
      P = Qxx + Kux.transpose()*Qux;
      
      dV[0] += ku.dot(Qu);
      dV[1] += ku.dot(Quu*ku/2);
    }
    
    //    if (debug)
    cout << "[I] Dmoc::Backward: current V=" << V << endl;    
  }
  
  
  template <typename T, int n, int c> 
    void Dmoc<T, n, c>::Forward() {

    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, n, n> Matrixnd;
    typedef Matrix<double, n, c> Matrixncd;
    typedef Matrix<double, c, n> Matrixcnd;
    typedef Matrix<double, c, c> Matrixcd;  
    
    // measured change in V
    double dVm = 1;
    
    double a = this->a;
    
    while (dVm > 0) {

      Vectornd dx = VectorXd::Zero(sys.n);
      T xn = xs[0];
      Vectorcd un;
      
      double Vm = 0;
      
      for (int k = 0; k < N; ++k) {
        const Vectorcd &u = us[k];
        Vectorcd &du = dus[k];
        
        const Vectorcd &ku = kus[k];
        const Matrixcnd &Kux = Kuxs[k];
        
        du = a*ku + Kux*dx;
        un = u + du;
        
        Rn<c> &U = (Rn<c>&)sys.U;
        if (U.bnd) {
          for (int j = 0; j < u.size(); ++j) 
            if (un[j] < U.lb[j]) {
              un[j] = U.lb[j];
              du[j] = un[j] - u[j];
            } else
              if (un[j] > U.ub[j]) {
                un[j] = U.ub[j];
                du[j] = un[j] - u[j];
              }
        }
        
        const double &t = ts[k];
        double h = ts[k+1] - t;

        double L = cost.L(t, xn, un, h);
        Vm += L;
        
        // update dx
        if (type == PURE) {
          dx = As[k]*dx + Bs[k]*du;
          cost.X.Retract(xn, xn, dx);
        } else {
          double h = ts[k+1] - t;
          T xn_(xn);
          sys.Step(xn_, t, xn, un, h);
          xn = xn_;
          cost.X.Lift(dx, xs[k+1], xn);

          //          cout << xn.gs[0] << " " << xn.r << " " << xn.vs[0] << " " << xn.dr << endl;

          //          cout << xs[k+1].gs[0] << " " << xs[k+1].r << " " << xs[k+1].vs[0] << " " << xs[k+1].dr << endl;
          assert(!std::isnan(dx[0]));
        }
      }
      
      double L = cost.L(ts[N], xn, un, 0);
      Vm += L;
      
      if (debug)
        cout << "[I] Dmoc::Forward: measured V=" << Vm << endl;
      
      dVm = Vm - V;
      
      if (dVm > 0) {
        a *= b1;
        if (a < 1e-12)
          break;
        if (debug)
          cout << "[I] Dmoc::Forward: step-size reduced a=" << a << endl;
        
        continue;
      }
      
      double r = dVm/(a*dV[0] + a*a*dV[1]);
      if (r < s1)
        a = b1*a;
      else 
        if (r >= s2) 
          a = b2*a;    
    
      if (debug)
        cout << "[I] Dmoc::Forward: step-size a=" << a << endl;    
    }
  }

  template <typename T, int n, int c> 
    void Dmoc<T, n, c>::Iterate() {
    
    Backward();
    Forward();
    for (int k = 0; k < N; ++k)
      us[k] += dus[k];
    Update();
  }
}

#endif
