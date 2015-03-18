#ifndef GCOP_PDDP_H
#define GCOP_PDDP_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <iostream>
#include "ddp.h"

namespace gcop {
  
  template <typename T, int n = Dynamic, int c = Dynamic, int np = Dynamic> 
    class PDdp : public Ddp<T, n, c, np> {
    
    typedef Eigen::Matrix<double, n, 1> Vectornd;
    typedef Eigen::Matrix<double, c, 1> Vectorcd;

    typedef Eigen::Matrix<double, n, n> Matrixnd;
    typedef Eigen::Matrix<double, n, c> Matrixncd;
    typedef Eigen::Matrix<double, c, n> Matrixcnd;
    typedef Eigen::Matrix<double, c, c> Matrixcd;  
    
    typedef Eigen::Matrix<double, np, 1> Vectormd;
    typedef Eigen::Matrix<double, np, np> Matrixmd;
    typedef Eigen::Matrix<double, n, np> Matrixnmd;
    typedef Eigen::Matrix<double, np, n> Matrixmnd;

    typedef Eigen::Matrix<double, c, np> Matrixcmd;
    typedef Eigen::Matrix<double, np, c> Matrixmcd;
    

  public:

    /**
     * Create a parameter-dependent optimal control problem using a system, a cost, and 
     * a trajectory given by a sequence of times, states, controls, and parameters.
     * The times ts must be given, the initial state xs[0] must be set, and
     * the controls us and parameters p will be used as an initial guess for the optimization.
     *
     * After initialization, every call to Iterate() will optimize the 
     * controls us, states xs, and parameters p and modify them accordingly. Problems involving
     * time-optimization will also modify the sequence of times ts.
     * 
     * @param sys system
     * @param cost cost
     * @param ts (N+1-vector) sequence of discrete times
     * @param xs (N+1-vector) sequence of discrete states
     * @param us (N-vector) sequence of control inputs
     * @param p (m-vector) parameters
     * @param dynParams how many of these are dynamic parameters (i.e. that the dynamics depends on), the first dynParams in the parameter vector p will be regarded as the dynamic parameters
     * @param update whether to update trajectory xs using initial state xs[0], inputs us, 
     *               and parameters p. This is necessary only if xs was not already generated.
     */    
    PDdp(System<T, n, c, np> &sys, Cost<T, n, c, np> &cost, 
          vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
          Vectormd &p, int dynParams = 0, bool update = true);
    
    /**
     * Perform one DDP iteration. Internally calls:
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


    int m;         ///< parameter space dimension
    int dynParams; ///< number of dynamic parameters (that the dynamics depend on)
    
    VectorXd &p;  ///< reference to parameters being optimized
    VectorXd dp;  ///< computed parameter increment      
  
    std::vector<Matrixnmd> Cs;

    std::vector<MatrixXd> Kuxs;

    VectorXd Lxa;
    MatrixXd Lxxa;
    Vectorcd Lu;
    Matrixcd Luu;
    MatrixXd Lxua;

    VectorXd v;
    MatrixXd P;
    
    double nu;    ///< current regularization factor nu
    double nu0;   ///< minimum regularization factor nu
    double dnu0;  ///< regularization factor modification step-size
    bool printDebug;


    bool pdX(const MatrixXd &P) {
      LLT<MatrixXd> llt;
      llt.compute(P);
      return llt.info() == Eigen::Success;
    }

  };
  
  
  using namespace std;
  using namespace Eigen;
  
  template <typename T, int n, int c, int np> 
    PDdp<T, n, c, np>::PDdp(System<T, n, c, np> &sys, 
                              Cost<T, n, c, np> &cost, 
                              vector<double> &ts, 
                              vector<T> &xs, 
                              vector<Matrix<double, c, 1> > &us,
                              Matrix<double, np, 1> &p, int dynParams, bool update) : 
    Ddp<T, n, c, np>(sys, cost, ts, xs, us, &p, false),
    m(p.size()), dynParams(dynParams), p(p), dp(m), Cs(this->N), 
    nu(1e-3), nu0(1e-3), dnu0(2) {
    
    Lxa.resize(sys.X.n + m);
    Lxxa.resize(sys.X.n + m, sys.X.n + m);
    Lxua.resize(sys.X.n + m, sys.U.n);
    Kuxs.resize(this->N);   
 
    for (int i = 0; i < this->N; ++i) {
      Kuxs[i] = MatrixXd::Zero(sys.U.n, sys.X.n + m);
      if (dynParams)
        Cs[i].resize(sys.X.n, dynParams);
    }


    if (update)
      Update();
    printDebug = false;
  }
  

  template <typename T, int n, int c, int np> 
    void PDdp<T, n, c, np>::Update(bool der) {
    
    VectorXd pdyn;
    if (dynParams)
      pdyn = this->p.head(dynParams);
    for (int k = 0; k < this->N; ++k) {
      double h = this->ts[k+1] - this->ts[k];
      if (der)
        this->sys.Step(this->xs[k+1], this->ts[k], this->xs[k], this->us[k], h, 
                       dynParams ? &pdyn : 0, 
                       &this->As[k], &this->Bs[k], dynParams ? &this->Cs[k] : 0);
      else
        this->sys.Step(this->xs[k+1], this->ts[k], this->xs[k], this->us[k], h, 
                       dynParams ? &pdyn : 0);
    }
  }
  
  template <typename T, int n, int c, int np> 
    void PDdp<T, n, c, np>::Backward() {
    
    int N = this->us.size();

    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, n, n> Matrixnd;
    typedef Matrix<double, n, c> Matrixncd;
    typedef Matrix<double, c, n> Matrixcnd; 
    typedef Matrix<double, c, c> Matrixcd;  
    typedef Matrix<double, Dynamic, 1> Vectormd;
    typedef Matrix<double, Dynamic, Dynamic> Matrixmd;
    typedef Matrix<double, n, Dynamic> Matrixnmd;
    typedef Matrix<double, Dynamic, n> Matrixmnd;
    
    double t = this->ts.back();
    const T &x = this->xs.back();
    const Vectorcd &u = this->us.back();


    Vectornd Lx;
    Matrixnd Lxx;
    Vectormd Lp;
    Matrixmd Lpp;
    Matrixmnd Lpx;
    Vectorcd Lu;
    Matrixcd Luu;
    Matrixncd Lxu;
   
    Lp.resize(m);
    Lpp.resize(m,m);
    Lpx.resize(m,n);

    double L = this->cost.L(t, x, u, 0, &p, 
                            &Lx, &Lxx, 0, 0, 0, 
                            &Lp, &Lpp, &Lpx);
    
    this->V = L;
    this->dV.setZero();

    /*
    cout << "Lxa: " << Lxa.rows() << "x" << Lxa.cols() << endl;
    cout << "Lx: " << Lx.rows() << "x" << Lx.cols() << endl;
    cout << "Lp: " << Lp.rows() << "x" << Lp.cols() << endl;
    cout << "Lxxa: " << Lxxa.rows() << "x" << Lxxa.cols() << endl;
    cout << "Lxx: " << Lxx.rows() << "x" << Lxx.cols() << endl;
    cout << "Lpp: " << Lpp.rows() << "x" << Lpp.cols() << endl;
    cout << "Lpx: " << Lpx.rows() << "x" << Lpx.cols() << endl;
    */

  
    Lxa.head(n) = Lx;
    Lxa.tail(m) = Lp;    
    Lxxa.block(0,0,n,n) = Lxx;
    Lxxa.block(n,n,m,m) = Lpp;
    Lxxa.block(n,0,m,n) = Lpx;
    Lxxa.block(0,n,n,m) = Lpx.transpose();

    /* 
    cout << "Lxa:" << endl << Lxa << endl;
    cout << "Lx:" << endl << Lx << endl;
    cout << "Lp:" << endl << Lp << endl;
    cout << "Lxxa:" << endl << Lxxa << endl;
    cout << "Lxx:" << endl << Lxx << endl;
    cout << "Lpp:" << endl << Lpp << endl;
    cout << "Lpx:" << endl << Lpx << endl;
    */

    v = Lxa;
    P = Lxxa;

    //cout << "v start: " << endl << v << endl;
    //cout << "P start: " << endl << P << endl;
    
    VectorXd Qx = VectorXd::Zero(this->sys.X.n + m);
    Vectorcd Qu;  
    
    MatrixXd Qxx = MatrixXd::Zero(this->sys.X.n + m, this->sys.X.n+m);
    Matrixcd Quu;
    Matrixcd Quum;
    MatrixXd Quxm =  MatrixXd::Zero(this->sys.U.n, this->sys.X.n+m);
    MatrixXd Qux = MatrixXd::Zero(this->sys.U.n, this->sys.X.n+m);

    MatrixXd Aat = MatrixXd::Zero(this->sys.X.n + m, this->sys.X.n+m);
    MatrixXd At = MatrixXd::Zero(this->sys.X.n, this->sys.X.n);
    MatrixXd Bt = MatrixXd::Zero(this->sys.U.n, this->sys.X.n);
    
    Matrixcd Ic = MatrixXd::Identity(this->sys.U.n, this->sys.U.n);

    for (int k = N - 1; k >=0; --k) {
      t = this->ts[k];
      const T &x = this->xs[k];
      const Vectorcd &u = this->us[k];
      double h = this->ts[k+1] - this->ts[k];
      assert(h >0);

      double L = this->cost.L(t, x, u, h, &p, &Lx, &Lxx, &Lu, &Luu, &Lxu, &Lp, &Lpp, &Lpx);
      Lxa.head(n) = Lx;
      Lxa.tail(m) = Lp;    
      Lxxa.block(0,0,n,n) = Lxx;
      Lxxa.block(n,n,m,m) = Lpp;
      Lxxa.block(n,0,m,n) = Lpx;
      Lxxa.block(0,n,n,m) = Lpx.transpose();
      Lxua.setZero();
      Lxua.block(0,0,n,c) = Lxu;

      //      cout << "L=" << L << endl;
      
      
      this->V += L;
      
      MatrixXd A = this->As[k];
      MatrixXd Aa = MatrixXd::Identity(this->sys.X.n + m, this->sys.X.n + m);
      Aa.block(0,0,n,n) = A;
      if(dynParams)
        Aa.block(0,n,n,dynParams) = this->Cs[k];

      MatrixXd B = this->Bs[k];//MatrixXd::Zero(this->sys.X.n + m, this->sys.U.n);
      //B.block(0,0,n,c) = this->Bs[k];

      Vectorcd &ku = this->kus[k];
      MatrixXd &Kux = Kuxs[k];
      
      Aat = Aa.transpose();
      At = A.transpose();
      Bt = B.transpose();
      
      Qx = Lxa + Aat*v;
      Qu = Lu + Bt*v.head(n);

      Qxx = Lxxa + Aat*P*Aa;
      Quu = Luu + Bt*P.block(0,0,n,n)*B;
      Qux.block(0,0,c,n) =  Bt*P.block(0,0,n,n)*A;     // assume Lux = 0
      Qux.block(0,n,c,m) =  Bt*P.block(0,n,n,m);     
      if(dynParams)
        Qux.block(0,n,c,dynParams) +=  Bt*P.block(0,0,n,n)*this->Cs[k];     

      double dmu = 1;
      
      LLT<Matrixcd> llt;
      
      printDebug = false;
     
      double mu = this->mu;
      if (this->debug) {
        if (!pdX(P)) {
        //  cout << "P[" << k << "] is not pd =" << endl << P << endl;
        //  cout << "Luu=" << endl << Luu << endl;
        }

        Matrixcd Pp = Bt*P.block(0,0,n,n)*B;
        if (!pdX(Pp)) {
        //  cout << "B'PB is not positive definite:" << endl << Pp << endl;
        }
      }

      if(false)
      {
        cout << "k=" << endl << k << endl;
        //cout << "x:" << endl << x << endl;
        cout << "u:" << endl << u << endl;
        cout << "p:" << endl << p << endl;
        cout << "Qxx:" << endl << Qxx << endl;
        cout << "Quu:" << endl << Quu << endl;
        cout << "Qux:" << endl << Qux << endl;
        cout << "P:" << endl << P << endl;
        cout << "v:" << endl << v << endl;
        cout << "Lxxa: " << endl << Lxxa << endl;
        cout << "Lxx: " << endl << Lxx << endl;
        cout << "Lpp: " << endl << Lpp << endl;
        cout << "Lpx: " << endl << Lpx << endl;
        cout << "Lxa: " << endl << Lxa << endl;
        cout << "Luu: " << endl << Luu << endl;
        cout << "Lu: " << endl << Lu << endl;
        cout << "Aa: " << endl << Aa << endl;
        cout << "B: " << endl << B << endl;
      }


      while (1) {
        Quum = Quu + mu*Ic;
        
        // From https://homes.cs.washington.edu/~todorov/papers/TassaIROS12.pdf
        //MatrixXd Pm = P+mu*MatrixXd::Identity(n+m,n+m);
        //Quum = Luu + Bt*Pm.block(0,0,n,n)*B;
        llt.compute(Quum);
        
        // if OK, then reduce mu and break
        if (llt.info() == Eigen::Success) {
          
          //Quxm.block(0,0,c,n) =  Bt*Pm.block(0,0,n,n)*A;     // assume Lux = 0
          //Quxm.block(0,n,c,m) =  Bt*Pm.block(0,n,n,m);     
          //if(dynParams)
          //  Quxm.block(0,n,c,dynParams) +=  Bt*Pm.block(0,0,n,n)*this->Cs[k];     

          // Tassa and Todorov recently proposed this quadratic rule, seems pretty good
          dmu = min(1/this->dmu0, dmu/this->dmu0);
          if (mu*dmu > this->mu0)
          { 
             mu = mu*dmu;
          }
          else
          {
             mu = this->mu0;
          }
          break;
        }
        else
        {
          // if negative then increase mu
          dmu = max(this->dmu0, dmu*this->dmu0);
          mu = max(this->mu0, mu*dmu);   
        }

        if (this->debug)
          cout << "[I] PDdp::Backward: increased mu=" << mu << " at k=" << k << endl;
        printDebug = true;

        if (mu > this->mumax) {
          cout << "[W] PDdp::Backward: mu= " << mu << " exceeded maximum !" << endl;
          if (this->debug)
            getchar();
          break;
        }

      }

      
      if(mu > this->mumax)
        break;


      ku = -llt.solve(Qu);
      Kux = -llt.solve(Qux);
      //Kux = -llt.solve(Quxm);

      assert(!std::isnan(ku[0]));
      assert(!std::isnan(Kux(0,0)));

      // Update used in https://homes.cs.washington.edu/~todorov/papers/TassaIROS12.pdf
      //v = Qx + Kux.transpose()*Quu*ku + Kux.transpose()*Qu + Qux.transpose()*ku;
      //P = Qxx + Kux.transpose()*Quu*Kux + Kux.transpose()*Qux + Qux.transpose()*Kux;
 
      v = Qx + Kux.transpose()*Qu;
      P = Qxx + Kux.transpose()*Qux;
      this->dV[0] += ku.dot(Qu);
      this->dV[1] += ku.dot(Quu*ku/2);

    }
    
    if (this->debug)
      cout << "[I] PDdp::Backward: current V=" << this->V << endl;
    
  }
  
  template <int n = Dynamic, int c = Dynamic> 
    Matrix<double, n, c> sym(const Matrix<double, n, c> &A) { return A+A.transpose(); };
  
  template <typename T, int n, int c, int np> 
    void PDdp<T, n, c, np>::Forward() {

    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, n, n> Matrixnd;
    typedef Matrix<double, n, c> Matrixncd;
    typedef Matrix<double, c, n> Matrixcnd;
    typedef Matrix<double, c, c> Matrixcd;  
    typedef Matrix<double, np, 1> Vectormd;
    typedef Matrix<double, np, np> Matrixmd;
    typedef Matrix<double, n, np> Matrixnmd;
    typedef Matrix<double, np, n> Matrixmnd;

    // find dp
    int N = this->us.size();

    LLT<Matrixmd> llt;
    double nu = this->nu;
    double dnu = 1;
    Vectormd bp = v.tail(m);

    while (1) {
      Matrixmd Apm = P.block(n,n,m,m);
      Apm.diagonal() = Apm.diagonal().array() + nu;
      
      llt.compute(Apm);
      // if OK, then reduce mu and break
      if (llt.info() == Eigen::Success) {
        dnu = min(1/this->dnu0, dnu/this->dnu0);
        if (nu*dnu > this->nu0)
          nu *= dnu;
        else
          nu = this->nu0;
        break;
      }
      
      // if negative then increase mu
      dnu = max(this->dnu0, dnu*this->dnu0);
      nu = max(this->nu0, nu*dnu);   
      if (this->debug)
      {
        cout << "[I] PDdp::Forward: increased nu=" << nu << endl;
      }
    }

    dp = -llt.solve(bp);

    // measured change in V
    double dVm = 1;      

    double a = this->a;

    Vectormd dp1 = dp;

    while (dVm > 0) {
      VectorXd dxa = VectorXd::Zero(this->sys.X.n + m);
      Vectornd dx = VectorXd::Zero(this->sys.X.n);

      dp = a*dp1;
      dxa.tail(m) = dp;
      //      Vectornd dx = Vectornd::Zero();
      T xn = this->xs[0];
      Vectorcd un;
      Vectormd pn = p + dp;
      
      VectorXd pdyn;
      if (dynParams)
        pdyn = pn.head(dynParams);
  
      double Vm = 0;
      
      for (int k = 0; k < N; ++k) {
        const Vectorcd &u = this->us[k];
        Vectorcd &du = this->dus[k];
        
        du = a*this->kus[k] + Kuxs[k]*dxa;
        un = u + du;
        
        //cout << "k=" << k << endl;
        //cout << "kus:" << endl << this->kus[k] << endl;
        //cout << "Kuxs:" << endl << Kuxs[k] << endl;

        const double &t = this->ts[k];
        double L = this->cost.L(t, xn, un, 0, &pn);
        Vm += L;
        
        // update dx
        if (this->type == this->PURE) {
          dx = this->As[k]*dx + this->Bs[k]*du;
          if(dynParams)
            dx += Cs[k]*dp.head(dynParams);
          this->sys.X.Retract(xn, xn, dx);
        } else {
          double h = this->ts[k+1] - t;
          T xn_(xn);
          this->sys.Step(xn_, t, xn, un, h, dynParams ? &pdyn : 0);
          xn = xn_;
          this->sys.X.Lift(dx, this->xs[k+1], xn);
        }

        dxa.head(n) = dx;
      }
      
      double L = this->cost.L(this->ts[N], xn, un, 0, &pn);
      Vm += L;
      
      if (this->debug)
        cout << "[I] PDdp::Forward: computed V=" << Vm << endl;
      
      dVm = Vm - this->V;
      
      if (dVm > 0) {
        a *= this->b1;
        if (a < 1e-12)
          break;
        if (this->debug)
          cout << "[I] PDdp::Forward: step-size reduced a=" << a << endl;
          //cout << "[I] PDdp::Forward: step-size reduced a=" << this->a << endl;
        
        
      }
      
      //double r = dVm/(this->a*this->dV[0] + this->a*this->a*this->dV[1]);
      double r = dVm/(a*this->dV[0] + a*a*this->dV[1]);
      if (r < this->s1)
        a *= this->b1;
      else 
        if (r >= this->s2)
          a *= this->b2;    
    
      if (this->debug)
        cout << "[I] PDdp::Forward: step-size a=" << a << endl;    
    }
  }

  template <typename T, int n, int c, int np> 
    void PDdp<T, n, c, np>::Iterate() {

    Backward();
    Forward();
    for (int k = 0; k < this->N; ++k)
      this->us[k] += this->dus[k];
    p = p + dp;
    Update();
  }
  
}

#endif

