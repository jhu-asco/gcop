// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_GDOCP_H
#define GCOP_GDOCP_H

#include "docp.h"
#include "ceres/ceres.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  /**
   * Gradient-based Discrete Optimal Control. 
   * Using the adjoint method for fast O(N) cost gradient computation. 
   * Uses the ceres library solver which supports 
   * STEEPEST_DESCENT, NONLINEAR_CONJUGATE_GRADIENT, BFGS and LBFGS
   * search directions and ARMIJO and WOLFE line search types.
   *
   * Author: Marin Kobilarov, marin(at)jhu.edu
   */
  template <typename T, int nx = Dynamic, int nu = Dynamic, int np = Dynamic> class GDocp :
    public Docp<T, nx, nu, np> {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, np, 1> Vectorpd;

    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
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
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    GDocp(System<T, nx, nu, np> &sys, Cost<T, nx, nu, np> &cost, 
          vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
          Vectorpd *p = 0, bool update = true);
    
    virtual ~GDocp();
    
    class GDocpCost : public ceres::FirstOrderFunction {
    public:
    GDocpCost(GDocp<T, nx, nu, np> &docp) : docp(docp) {};

      virtual ~GDocpCost() {};

      virtual bool Evaluate(const double* parameters,
                            double* cost,
                            double* gradient) const {
        assert(parameters);
        for (int i = 0; i < docp.us.size(); ++i)
          memcpy(docp.us[i].data(), parameters + i*docp.sys.U.n, docp.sys.U.n*sizeof(double));
        docp.Update();   // forward
        if (gradient) {
          if (docp.prec)
            docp.EvaluatePrec(); // backward
          else
            docp.Evaluate(); // backward            
        } else {
          docp.J = docp.ComputeCost(); // backward          
        } 
      
        if (cost)
          *cost = docp.J;
        if (gradient)
          for (int i = 0; i < docp.us.size(); ++i)
            memcpy(gradient + i*docp.sys.U.n, docp.Jus[i].data(), docp.sys.U.n*sizeof(double));
        return true;
      }

      virtual int NumParameters() const { return docp.us.size()*docp.sys.U.n; }

      GDocp<T, nx, nu, np> &docp;
    };


    /**
     * Perform one GDOCP iteration. Internally calls:
     *
     */
    void Iterate();

    /**
     *  Evaluate cost and its gradient
     */
    void Evaluate();
    

    /**
     *  Evaluate cost and its gradient, preconditioned using the
     * inverse of the Hessian matrix
     */
    void EvaluatePrec();

    int N;        ///< number of discrete trajectory segments    
    
    double mu;    ///< current regularization factor mu
    double mu0;   ///< minimum regularization factor mu
    double dmu0;  ///< regularization factor modification step-size
    double mumax; ///< maximum regularization factor mu
    double a;     ///< step-size

    Vectornd Lx;
    Matrixnd Lxx;
    Vectorcd Lu;
    Matrixcd Luu;

    Matrixnd P;
    Vectornd v;          ///< multiplier

    vector<Vectorcd> Jus; ///< gradient of J with respect to controls

    bool prec;    ///< use Hessian-approximation preconditioning
    
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
    GDocp<T, nx, nu, np>::GDocp(System<T, nx, nu, np> &sys, 
                            Cost<T, nx, nu, np> &cost, 
                            vector<double> &ts, 
                            vector<T> &xs, 
                            vector<Matrix<double, nu, 1> > &us,
                            Matrix<double, np, 1> *p,
                            bool update) : Docp<T, nx, nu, np>(sys, cost, ts, xs, us, p, update),
    N(us.size()), 
    mu(1e-3), mu0(1e-3), dmu0(2), mumax(1e6), a(1), 
    Jus(us.size()),
    prec(false)
    {
      assert(N > 0);
      assert(sys.X.n > 0);
      assert(sys.U.n > 0);

      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      assert(us.size() == N);

      if (nu == Dynamic) {
        for (int i = 0; i < N; ++i) {
          Jus[i].resize(sys.U.n);
        }
        Lu.resize(sys.U.n);
        Luu.resize(sys.U.n, sys.U.n);
      }
      
      if (nx == Dynamic) {
        Lx.resize(sys.X.n);
        Lxx.resize(sys.X.n, sys.X.n);
        v.resize(sys.X.n);       
        P.resize(sys.X.n, sys.X.n);
      }

      if (update) {
        this->Update();
      }
    }
  
  template <typename T, int nx, int nu, int np> 
    GDocp<T, nx, nu, np>::~GDocp()
    {
    }


  template <typename T, int nx, int nu, int np> 
    void GDocp<T, nx, nu, np>::EvaluatePrec() {

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
    
    this->J = L;
    
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
      
      this->J += L;
      
      const Matrixnd &A = this->As[k];
      const Matrixncd &B = this->Bs[k];
      
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
      
      Vectorcd ku = -llt.solve(Qu);
      Matrixcnd Kux = -llt.solve(Qux);
      
      assert(!std::isnan(ku[0]));
      assert(!std::isnan(Kux(0,0)));

      v = Qx + Kux.transpose()*Qu;
      P = Qxx + Kux.transpose()*Qux;
      
      Jus[k] = -ku;
      //      dV[0] += ku.dot(Qu);
      //      dV[1] += ku.dot(Quu*ku/2);
    }
    
    //    if (debug)
    cout << "[I] Ddp::Backward: current J=" << this->J << endl;    
  }
    
  template <typename T, int nx, int nu, int np> 
    void GDocp<T, nx, nu, np>::Evaluate() {

    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    
    double t = this->ts.back();
    const T &x = this->xs.back();
    const Vectorcd &u = this->us.back();
    double L = this->cost.L(t, x, u, 0, 0, &Lx, 0);
    
    this->J = L;    
    v = Lx;
    
    for (int k = N-1; k >=0; --k) {
      
      t = this->ts[k];
      double h = this->ts[k+1] - t;

      const T &x = this->xs[k];
      const Vectorcd &u = this->us[k];
      double L = this->cost.L(t, x, u, h, 0, &Lx, 0, &Lu, 0);

      if (std::isnan(L)) {
        cout << "NAN" << " k=" << k << " Lx=" << Lx << endl;
        getchar();
      }
      
      this->J += L;
      Jus[k] = Lu + this->Bs[k].transpose()*v;
      v = Lx + this->As[k].transpose()*v;
    }    
    //    if (debug)
    cout << "[I] GDocp:: current J=" << this->J << endl;    
  }

  template <typename T, int nx, int nu, int np> 
    void GDocp<T, nx, nu, np>::Iterate() {
    
    double parameters[this->us.size()*this->sys.U.n];
    for (int i = 0; i < this->us.size(); ++i)
      memcpy(parameters + i*this->sys.U.n, this->us[i].data(), this->sys.U.n*sizeof(double));
    
    
    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout = true;
    //    options.line_search_direction_type = ceres::NONLINEAR_CONJUGATE_GRADIENT;

    //    options.line_search_type = ceres::ARMIJO;
    //    options.line_search_direction_type = ceres::STEEPEST_DESCENT;

    options.max_num_iterations = 200;
    ceres::GradientProblemSolver::Summary summary;
    ceres::GradientProblem problem(new GDocpCost(*this));
    ceres::Solve(options, problem, parameters, &summary);        
    cout << summary.FullReport();

    for (int i = 0; i < this->us.size(); ++i)
      memcpy(this->us[i].data(), parameters + i*this->sys.U.n, this->sys.U.n*sizeof(double));
    this->Update(false);    
  }
}

#endif
