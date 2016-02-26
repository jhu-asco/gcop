// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_DDP_H
#define GCOP_DDP_H

#include "docp.h"
#include <exception>      // std::exception
#include <stdexcept>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  /**
   * Differential Dynamic Programming
   *
   * Author: Marin Kobilarov 
   */
  template <typename T, int nx = Dynamic, int nu = Dynamic, int np = Dynamic> class Ddp :
    public Docp<T, nx, nu, np> {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, np, 1> Vectorpd;

    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;
    
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
     * @param p parameters (optional)
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    Ddp(System<T, nx, nu, np> &sys, Cost<T, nx, nu, np> &cost, 
        vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
        Vectorpd *p = 0, bool update = true);


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
     * @param traj trajectory
     * @param cost cost
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    Ddp(Trajectory<T, nx, nu, np> &traj, Cost<T, nx, nu, np> &cost, bool update = true);
    
    virtual ~Ddp();
   
    /**
     * Perform one DDP iteration. Internally calls:
     *
     * Backward -> Forward -> Update. The controls us and trajectory xs
     * are updated. 
     */
    void Iterate();

    /**
     *  Forward pass
     */
    void Forward();

    /**
     *  Backward pass
     */
    void Backward();
        
    int N;        ///< number of discrete trajectory segments
    
    double mu;    ///< current regularization factor mu
    double mu0;   ///< minimum regularization factor mu
    double dmu0;  ///< regularization factor modification step-size
    double mumax; ///< maximum regularization factor mu
    double a;     ///< step-size
    
    std::vector<Vectorcd> dus;  ///< computed control change
    
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

    double sigma; ///< Armijo sigma factor (0.1 by default)
    double beta;  ///< Armijo beta factor  (0.25 by default)

    double amin;  ///< minimum step-size: end line-search if a<amin (default is 1e-10)
    double cVmin; ///< cost-change minimum: end line-search if V_a - V < cVmin (default is 1-10)
    
    double s1;   ///< Armijo-Goldstein step-size control factor s1
    double s2;   ///< Armijo-Goldstein step-size control factor s2
    double b1;   ///< Armijo-Goldstein step-size control factor b1
    double b2;   ///< Armijo-Goldstein step-size control factor b2
    
    char type;   ///< type of algorithm (choices are PURE, DDP, LQS), LQS is default. In the current implementation second-order expansion of the dynamics is ignored, so DDP and LQS are identical. Both LQS and DDP are implemented using the true nonlinear dynamics integration in the Forward step, while PURE uses the linearized version in order to match exactly the true Newton step. 
    
    static const char PURE = 0;     ///< PURE version of algorithm (i.e. stage-wise Newton)
    static const char DDP = 1;      ///< DDP version of algorithm
    static const char LQS = 2;      ///< Linear-Quadratic Subproblem version of algorithm due to Dreyfus / Dunn&Bertsekas / Pantoja
    
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
  
  template <typename T, int nx, int nu, int np> 
    Ddp<T, nx, nu, np>::Ddp(Trajectory<T, nx, nu, np> &traj, 
                            Cost<T, nx, nu, np> &cost, 
                            bool update) : Ddp(traj.sys, cost, traj.ts, traj.xs, traj.us, traj.p, update) {
  }

  template <typename T, int nx, int nu, int np> 
    Ddp<T, nx, nu, np>::Ddp(System<T, nx, nu, np> &sys, 
                            Cost<T, nx, nu, np> &cost, 
                            vector<double> &ts, 
                            vector<T> &xs, 
                            vector<Matrix<double, nu, 1> > &us,
                            Matrix<double, np, 1> *p,
                            bool update) : Docp<T, nx, nu, np>(sys, cost, ts, xs, us, p, update),
    N(us.size()), 
    mu(1e-3), mu0(1e-3), dmu0(2), mumax(1e6), a(1), 
    dus(N), kus(N), Kuxs(N), 
    sigma(0.1), beta(.25), amin(1e-10), cVmin(1e-10),
    s1(0.1), s2(0.5), b1(0.25), b2(2),
    type(LQS)
    {
      assert(N > 0);
      assert(sys.X.n > 0);
      assert(sys.U.n > 0);

      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      assert(us.size() == N);

      if (nx == Dynamic || nu == Dynamic) {
        for (int i = 0; i < N; ++i) {
          dus[i].resize(sys.U.n);
          //          As[i].resize(sys.X.n, sys.X.n);
          //          Bs[i].resize(sys.X.n, sys.U.n);
          kus[i].resize(sys.U.n);
          Kuxs[i].resize(sys.U.n, sys.X.n);
        }

        Lx.resize(sys.X.n);
        Lxx.resize(sys.X.n, sys.X.n);
        Lu.resize(sys.U.n);
        Luu.resize(sys.U.n, sys.U.n);

        P.resize(sys.X.n, sys.X.n);
        v.resize(sys.X.n);       
      }

      if (update) {
        this->Update();
      }
    }
  
  template <typename T, int nx, int nu, int np> 
    Ddp<T, nx, nu, np>::~Ddp()
    {
    }
    
  template <typename T, int nx, int nu, int np> 
    void Ddp<T, nx, nu, np>::Backward() {

    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  
    
    double t = this->ts.back();
    const T &x = this->xs.back();
    const Vectorcd &u = this->us.back();
    double L = this->cost.L(t, x, u, 0, 0, &Lx, &Lxx);
    
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
    
    Matrixcd Ic = MatrixXd::Identity(this->sys.U.n, this->sys.U.n);
    
    for (int k = N-1; k >=0; --k) {
      
      t = this->ts[k];
      double h = this->ts[k+1] - t;

      const T &x = this->xs[k];
      const Vectorcd &u = this->us[k];
      double L = this->cost.L(t, x, u, h, 0, &Lx, &Lxx, &Lu, &Luu);

      if (std::isnan(L)) {
        cout << "NAN" << " k=" << k << " Lx=" << Lx << " Lxx=" << Lxx << endl;
        getchar();
      }
      
      V += L;
      
      const Matrixnd &A = this->As[k];
      const Matrixncd &B = this->Bs[k];
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
     
      if (this->debug) {
      if (!pd(P)) {
        cout << "P[" << k << "] is not pd =" << endl << P << endl;
        cout << "Luu=" << endl << Luu << endl;
      }

      Matrixcd Pp = Bt*P*B;
      if (!pdX(Pp)) {
        cout << "B'PB is not positive definite:" << endl << Pp << endl;
      }
      }

      LLT<Matrixcd> llt;     
      
      while (1) {
        Quum = Quu + mu*Ic;
        llt.compute(Quum);
        
        // if OK, then reduce mu and break
        if (llt.info() == Eigen::Success) {
          // this is the standard quadratic rule specified by Tassa and Todorov
          dmu = min(1/this->dmu0, dmu/this->dmu0);
          if (mu*dmu > this->mu0)
            mu = mu*dmu;
          else
            mu = this->mu0;
          
        if (this->debug) 
          cout << "[I] Ddp::Backward: reduced mu=" << mu << " at k=" << k << endl;


          break;
        }
        
        // if negative then increase mu
        dmu = max(this->dmu0, dmu*this->dmu0);
        mu = max(this->mu0, mu*dmu);
        
        if (this->debug) {
          cout << "[I] Ddp::Backward: increased mu=" << mu << " at k=" << k << endl;
          cout << "[I] Ddp::Backward: Quu=" << Quu << endl;
        }

        if (mu > mumax) {
          cout << "[W] Ddp::Backward: mu= " << mu << " exceeded maximum !" << endl;          
          if (this->debug)
            getchar();
          break;
        }
      }

      if (mu > mumax)
        break;
      
      ku = -llt.solve(Qu);
      Kux = -llt.solve(Qux);
      
      assert(!std::isnan(ku[0]));
      assert(!std::isnan(Kux(0,0)));

      v = Qx + Kux.transpose()*Qu;
      P = Qxx + Kux.transpose()*Qux;
      
      dV[0] += ku.dot(Qu);
      dV[1] += ku.dot(Quu*ku/2);
    }
    
    if (this->debug)
    cout << "[I] Ddp::Backward: current V=" << V << endl;    
  }

  template <typename T, int nx, int nu, int np> 
    void Ddp<T, nx, nu, np>::Forward() {

    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  
    
    // measured change in V
    double cV = 1;
    
    double a = this->a;

    bool acc = false;
    
    while (1) {

      Vectornd dx = VectorXd::Zero(this->sys.X.n);
      dx.setZero();//Redundancy
      T xn = this->xs[0];
      this->sys.Reset(xn,this->ts[0]);//Reset to initial state
      Vectorcd un;
      
      double Vm = 0;
      
      for (int k = 0; k < N; ++k) {
        const Vectorcd &u = this->us[k];
        Vectorcd &du = dus[k];
        
        const Vectorcd &ku = kus[k];
        const Matrixcnd &Kux = Kuxs[k];
        
        du = a*ku + Kux*dx;
        un = u + du;
				//[DEBUG]:
        //cout<<"du[ "<<k<<"]:\t"<<du.transpose()<<endl;
	//			cout<<"dx[ "<<k<<"]:\t"<<dx.transpose()<<endl;
	//			cout<<"ku[ "<<k<<"]:\t"<<ku.transpose()<<endl;
	//			cout<<"Kux[ "<<k<<"]:"<<endl<<Kux<<endl;
        
        Rn<nu> &U = (Rn<nu>&)this->sys.U;
        if (U.bnd) {
          U.Bound(un, du, u);
        }
        
        const double &t = this->ts[k];
        double h = this->ts[k+1] - t;

        double L = this->cost.L(t, xn, un, h);
        Vm += L;
        
        // update dx
        if (type == PURE) {
          dx = this->As[k]*dx + this->Bs[k]*du;
          this->sys.X.Retract(xn, xn, dx);
        } else {
          double h = this->ts[k+1] - t;
          //Adding nan catching :
          try {
            this->sys.Step(xn, un, h, this->p);
          }
          catch(std::exception &e)
            {
              //[DEBUG] Statement
              std::cerr << "exception caught: " << e.what() << '\n';
              std::cout<<" Iteration counter: "<<k<<endl;
              std::cout<<" u: "<<u.transpose()<<endl;
              std::cout<<" du: "<<du.transpose()<<endl;
              //More debug statements:
              //std::cout<<" a: "<<a<<"\t ku: "<<ku.transpose()<<"\t dx: "<<dx.transpose()<<"\n kux: \n"<<Kux<<endl;
              //[Culprit dx]
              
              throw std::runtime_error(std::string("Nan observed"));
              return;
            }
          //cout<<"dx: "<<dx.transpose()<<endl;//[DEBUG]//Previous iteration dx
          //std::cout<<" du: "<<du.transpose()<<endl;//[DEBUG]
          
          this->sys.X.Lift(dx, this->xs[k+1], xn);
          
          //          cout<<"xs[ "<<(k+1)<<"]:\t"<< ((Body2dState&)this->xs[k+1]).first<< " " <<((Body2dState&)this->xs[k+1]).second << endl;

          //          cout<<"un[ "<<(k+1)<<"]:\t"<<un.transpose()<<endl;
          
          //          cout << xn.gs[0] << " " << xn.r << " " << xn.vs[0] << " " << xn.dr << endl;
          
          //          cout << this->xs[k+1].gs[0] << " " << this->xs[k+1].r << " " << this->xs[k+1].vs[0] << " " << this->xs[k+1].dr << endl;
          assert(!std::isnan(dx[0]));
        }
      }
      
      double L = this->cost.L(this->ts[N], xn, un, 0);
      Vm += L;
      cV = Vm - V;
      
      if (this->debug)
        cout << "[I] Ddp::Forward: measured V=" << Vm << endl;
      
      //      cout << "V-Vm=" << V-Vm << " dV[0]=" << dV[0] << " a=" << a << endl;
     
      if (a < amin || fabs(cV) < cVmin)
        break;

      if (V - Vm < -sigma*a*dV[0]) {
        a *= beta;
        continue;
      } else {
        break;
      }

      
      // the rest below is currently not activated

      // if step is still not acceptable
      if (acc)
        break;            

      // Armijo-Goldstein
      if (cV > 0) {
        a *= b1;
        if (a < 1e-12)
          break;
        if (this->debug)
          cout << "[I] Ddp::Forward: step-size reduced a=" << a << endl;

        acc = false;
        continue;
      }

      double r = cV/(a*dV[0] + a*a*dV[1]);
      //      double r = cV/(a*dV[0]);
      if (r < s1)
        a = b1*a;
      else 
        if (r >= s2) 
          a = b2*a;
        else
          break;

      if (a < 1e-12)
        break;
      
      if (this->debug)
        cout << "[I] Ddp::Forward: step-size a=" << a << endl;    
    }
    this->J = V + cV;//Set the optimal cost after one iteration
  }

  template <typename T, int nx, int nu, int np> 
    void Ddp<T, nx, nu, np>::Iterate() {
    
    Backward();
    Forward();
    for (int k = 0; k < N; ++k)
      this->us[k] += dus[k];
    this->Update();
  }
}

#endif
