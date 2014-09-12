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

    std::vector<Vectornd> Lxs;
    std::vector<Matrixnd> Lxxs;
    std::vector<Vectorcd> Lus;
    std::vector<Matrixcd> Luus;
    std::vector<Matrixncd> Lxus;
    std::vector<Vectormd> Lps;
    std::vector<Matrixmd> Lpps;
    std::vector<Matrixmnd> Lpxs;
    
    std::vector<Matrixcmd> Kups;

    std::vector<Matrixnmd> Hxs;
    std::vector<Vectornd> hxs;
    
    Vectormd Lp;
    Matrixmd Lpp;
    Matrixmnd Lpx;    

    double nu;    ///< current regularization factor nu
    double nu0;   ///< minimum regularization factor nu
    double dnu0;  ///< regularization factor modification step-size

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
    Lxs(this->N+1), Lxxs(this->N+1), 
    Lus(this->N), Luus(this->N), Lxus(this->N),    
    Lps(this->N+1), Lpps(this->N+1), Lpxs(this->N+1),
    Kups(this->N),
    Hxs(this->N+1), hxs(this->N+1),
    nu(1e-3), nu0(1e-3), dnu0(2) {
    
    for (int k = 0; k <= this->N; ++k) {
      if (k < this->N) {
        if (dynParams)
          Cs[k].resize(sys.X.n, dynParams);
        Kups[k].resize(sys.U.n, m);
      }
      Hxs[k].resize(sys.X.n, m);
      
      Lps[k].resize(m, 1);
      Lpps[k].resize(m, m);
      Lpxs[k].resize(m, sys.X.n);
    }


    if (n == Dynamic || c == Dynamic) {      
      for (int i = 0; i <= this->N; ++i) {
        Lxs[i].resize(sys.X.n);
        Lxxs[i].resize(sys.X.n, sys.X.n);

        if (i < this->N) {
          Lus[i].resize(sys.U.n);
          Luus[i].resize(sys.U.n, sys.U.n);
          Lxus[i].resize(sys.X.n, sys.U.n);
        }

        hxs[i].resize(sys.X.n);
      }
    }

    if (update)
      Update();
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

    double L = this->cost.L(t, x, u, 0, &p, 
                            &Lxs[N], &Lxxs[N], 0, 0, 0, 
                            &Lps[N], &Lpps[N], &Lpxs[N]);
    
    this->V = L;
    this->dV.setZero();
    
    Vectornd v = Lxs[N];
    Matrixnd P = Lxxs[N];
    Matrixmnd D = Lpxs[N];
    
    Vectornd Qx;
    Vectorcd Qu;  
    Vectormd Qp;
    
    Matrixnd Qxx;
    Matrixcd Quu;
    Matrixcd Quum;
    Matrixcnd Qux;
    Matrixmcd Qpu;

    Matrixmnd Qpx;
    
    Matrixnd At;
    Matrixcnd Bt;
    Matrixmnd Ct;
    
    Matrixcd Ic = MatrixXd::Identity(this->sys.U.n, this->sys.U.n);

    for (int k = N - 1; k >=0; --k) {
      t = this->ts[k];
      const T &x = this->xs[k];
      const Vectorcd &u = this->us[k];
      
      Vectornd  &Lx = Lxs[k];
      Matrixnd &Lxx = Lxxs[k];
      Vectorcd  &Lu = Lus[k];
      Matrixcd &Luu = Luus[k];
      Matrixncd &Lxu = Lxus[k];
      Vectormd &Lp = Lps[k];
      Matrixmd &Lpp = Lpps[k];
      Matrixmnd &Lpx = Lpxs[k];

      double h = this->ts[k+1] - this->ts[k];
      assert(h >0);
      double L = this->cost.L(t, x, u, h, &p, &Lx, &Lxx, &Lu, &Luu, &Lxu, &Lp, &Lpp, &Lpx);

      //      cout << "L=" << L << endl;
      
      /*
      cout << "Lx=" << Lx << endl;
      cout << "Lxx=" << Lxx << endl;
      cout << "Lu=" << Lu << endl;
      cout << "Luu=" << Luu << endl;
      cout << "Lxu=" << Lxu << endl;
      cout << "Lp=" << Lp << endl;
      cout << "Lpp=" << Lpp << endl;
      cout << "Lpx=" << Lpx << endl;
      */
      
      this->V += L;
      
      const Matrixnd &A = this->As[k];
      const Matrixncd &B = this->Bs[k];
      const Matrixnmd &C = Cs[k];

      Vectorcd &ku = this->kus[k];
      Matrixcnd &Kux = this->Kuxs[k];
      Matrixcmd &Kup = Kups[k];
      
      At = A.transpose();
      Bt = B.transpose();
      Ct = C.transpose();
      
      Qx = Lx + At*v;
      Qu = Lu + Bt*v;
      Qp = Lp;
      if (dynParams)
        Qp.head(dynParams) += Ct*v;

      Qxx = Lxx + At*P*A;
      Quu = Luu + Bt*P*B;
      Qux = Bt*P*A;     // assume Lux = 0
      Qpx = Lpx + D*A;
      if (dynParams)
        Qpx.topRows(dynParams) += Ct*P*A;

      Qpu = D*B;  //assume Lpu=0
      if (dynParams)
        Qpu.topRows(dynParams) += Ct*P*B;

      /*
      Qpx = Lpx;
      if (dynParams)
        Qpx.topRows<dynParams>() += Ct*P*A;
      */

      double mu = this->mu;
      double dmu = 1;
      
      LLT<Matrixcd> llt;
      
      while (1) {
        Quum = Quu + this->mu*Ic;
        llt.compute(Quum);
        
        // if OK, then reduce mu and break
        if (llt.info() == Eigen::Success) {
          // Tassa and Todorov recently proposed this quadratic rule, seems pretty good
          dmu = min(1/this->dmu0, dmu/this->dmu0);
          if (this->mu*dmu > this->mu0)
            this->mu *= dmu;
          else
            this->mu = this->mu0;
          break;
        }
        
        // if negative then increase mu
        dmu = max(this->dmu0, dmu*this->dmu0);
        this->mu = max(this->mu0, this->mu*dmu);   
        
        if (this->debug)
          cout << "[I] PDdp::Backward: increased mu=" << mu << " at k=" << k << endl;
      }
      
      ku = -llt.solve(Qu);
      Kux = -llt.solve(Qux);
      Kup = -llt.solve(Qpu.transpose());

      
      v = Qx + Kux.transpose()*Qu;//% + Kup.transpose()*Qp;
      P = Qxx + Kux.transpose()*Qux;// + Kup.transpose()*Qpx;
      
      D = Qpx + Qpu*Kux;
      //      dV[0] += ku.dot(Qu);
      //      dV[1] += ku.dot(Quu*ku/2);
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

    hxs[0].setZero();
    Hxs[0].setZero();

    Matrixmd Ap = Lpps[0] + Kups[0].transpose()*this->Luus[0]*Kups[0];
    Vectormd bp = Lps[0] + Kups[0].transpose()*(this->Lus[0] + this->Luus[0]*this->kus[0]);

    Matrixnd E;
    for (int k = 1; k <= N; ++k) {
      E = this->As[k-1] + this->Bs[k-1]*this->Kuxs[k-1];

      Hxs[k] = E*Hxs[k-1] + this->Bs[k-1]*Kups[k-1];
      if (dynParams)
        Hxs[k].leftCols(dynParams) += Cs[k-1];
      hxs[k] = E*hxs[k-1] + this->Bs[k-1]*this->kus[k-1];

      // assume Lup = 0

      Ap = Ap + Hxs[k].transpose()*this->Lxxs[k]*Hxs[k] + Lpxs[k]*Hxs[k] + (Lpxs[k]*Hxs[k]).transpose() + Lpps[k]; 

      if (k < N)
        Ap = Ap + (this->Kuxs[k]*Hxs[k] + Kups[k]).transpose()*this->Luus[k]*(this->Kuxs[k]*Hxs[k] + Kups[k]);
      
 
      bp = bp + Hxs[k].transpose()*(this->Lxs[k] + this->Lxxs[k]*hxs[k]) + Lpxs[k]*hxs[k] + Lps[k]; 

      if (k<N)
        bp = bp + (this->Kuxs[k]*Hxs[k] + Kups[k]).transpose()*(this->Lus[k] + this->Luus[k]*(this->Kuxs[k]*hxs[k] + this->kus[k])); 
 
    }

    //    cout << "Ap=" << Ap << endl;
    //    cout << "bp=" << bp << endl;

    //    bool pd = true;

    //    if (pd) {
    
    LLT<Matrixmd> llt;
    double nu = this->nu;
    double dnu = 1;
    
    while (1) {
      Matrixmd Apm = Ap;
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
        cout << "[I] PDdp::Backward: increased nu=" << nu << endl;
    }
    
    dp = -llt.solve(bp);

    //    dp = -dp./diag(Ap);

    //      cout << "[W] Pddp::Backward: Ap not positive definite!" << endl;
    //      ColPivHouseholderQR<Matrixmd> dec(Ap);
    //      dp = -dec.solve(bp);    

    this->dV.setZero();

    // compute expected cost
    for (int k = 0; k <= N; ++k) {
      Vectornd dx = Hxs[k]*dp + hxs[k];
      
      this->dV[0] += Lps[k].dot(dp) + Lxs[k].dot(dx);
      this->dV[1] += dp.dot(Lpps[k]*dp)/2 + dp.dot(Lpxs[k]*dx) + dx.dot(Lxxs[k]*dx)/2;
      
      if (k < N) {
        Vectorcd du = this->Kuxs[k]*dx + Kups[k]*dp + this->kus[k];
        this->dV[0] += Lus[k].dot(du);
        this->dV[1] += du.dot(Luus[k]*du)/2;
      }
    }

    // measured change in V
    double dVm = 1;
    
    double a = this->a;

    Vectormd dp1 = dp;

    while (dVm > 0) {
      Vectornd dx = VectorXd::Zero(this->sys.X.n);
			dx.setZero();//Redundancy
      //      Vectornd dx = Vectornd::Zero();
      T xn = this->xs[0];
      Vectorcd un;
      
      double Vc = 0;
      
      dp = a*dp1;
      Vectormd pn = p + dp;
      VectorXd pdyn;
      if (dynParams)
        pdyn = pn.head(dynParams);

      for (int k = 0; k < N; ++k) {
        const Vectorcd &u = this->us[k];
        Vectorcd &du = this->dus[k];
        
        du = a*this->kus[k] + this->Kuxs[k]*dx + Kups[k]*dp;
        un = u + du;
        
        const double &t = this->ts[k];
        double L = this->cost.L(t, xn, un, 0, &pn);
        Vc += L;
        
        // update dx
        if (this->type == this->PURE) {
          dx = this->As[k]*dx + this->Bs[k]*du;
          if (dynParams)
            dx += Cs[k]*dp.head(dynParams);
          this->sys.X.Retract(xn, xn, dx);
        } else {
          double h = this->ts[k+1] - t;
          T xn_(xn);
          this->sys.Step(xn_, t, xn, un, h, dynParams ? &pdyn : 0);
          xn = xn_;
          this->sys.X.Lift(dx, this->xs[k+1], xn);
        }
      }
      
      double L = this->cost.L(this->ts[N], xn, un, 0, &pn);
      Vc += L;
      
      if (this->debug)
        cout << "[I] PDdp::Forward: computed V=" << Vc << endl;
      
      dVm = Vc - this->V;
      
      if (dVm > 0) {
        a *= this->b1;
        if (a < 1e-12)
          break;
        if (this->debug)
          cout << "[I] PDdp::Forward: step-size reduced a=" << this->a << endl;
        
        continue;
      }
      
      double r = dVm/(this->a*this->dV[0] + this->a*this->a*this->dV[1]);
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

