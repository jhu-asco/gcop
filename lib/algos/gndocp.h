#ifndef GCOP_GNDOCP_H
#define GCOP_GNDOCP_H

#include "docp.h"
#include "lqcost.h"
#include "tparam.h"

#include <unsupported/Eigen/NonLinearOptimization>

#ifdef GCOP_GNDOCP_CERES
#include "ceres/ceres.h"
#include "glog/logging.h"
#endif

//#define USE_SAMPLE_NUMERICAL_DIFF
/*#ifdef USE_SAMPLE_NUMERICAL_DIFF
#include "samplenumericaldiff.h"
#endif
*/

namespace gcop {

#ifdef GCOP_GNDOCP_CERES
using ceres::DynamicNumericDiffCostFunction;
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
#endif
  
  using namespace std;
  using namespace Eigen;

 
  template <typename T, 
    int _nx, 
    int _nu,
    int _np,
    int _ng,
    int _ntp> class GnCost;
  
  template <typename T, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic,
    int _ng = Dynamic,
    int _ntp = Dynamic> class GnDocp : 
    public Docp<T, _nx, _nu, _np> {
    
    typedef Matrix<double, _ng, 1> Vectorgd;
    typedef Matrix<double, _ng, _nx> Matrixgxd;
    typedef Matrix<double, _ng, _nu> Matrixgud;
    typedef Matrix<double, _ng, _np> Matrixgpd;
    
    typedef Matrix<double, _nx, 1> Vectornd;
    typedef Matrix<double, _nu, 1> Vectorcd;

    typedef Matrix<double, _nx, _nx> Matrixnd;
    typedef Matrix<double, _nx, _nu> Matrixncd;
    typedef Matrix<double, _nu, _nx> Matrixcnd;
    typedef Matrix<double, _nu, _nu> Matrixcd;
    
    typedef Matrix<double, _np, 1> Vectormd;
    
    typedef Matrix<double, _np, _np> Matrixmd;
    typedef Matrix<double, _nx, _np> Matrixnmd;
    typedef Matrix<double, _np, _nx> Matrixmnd;    
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
     * @param tparam trajectory parametrization
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param p static parameter vector
     * @param ng 
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    GnDocp(System<T, _nx, _nu, _np> &sys, 
           LqCost<T, _nx, _nu, _np, _ng> &cost, Tparam<T, _nx, _nu, _np, _ntp> &tparam,
           vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
           Vectormd *p = 0, bool update = true);
    
    virtual ~GnDocp();

    /**
     * Perform one DOCP iteration. Internally calls:
     * are updated. 
     */
    void Iterate();

    void Reset();
    Tparam<T, _nx, _nu, _np, _ntp> &tparam;

    int info;
    double fnorm, covfac;

    int inputs;
    int values;

    double numdiff_stepsize;   ///< The step size for perturbations
     
    VectorXd s;  ///< optimization vector 
    
    GnCost<T, _nx, _nu, _np, _ng, _ntp> *functor;
//#ifndef USE_SAMPLE_NUMERICAL_DIFF
    NumericalDiff<GnCost<T, _nx, _nu, _np, _ng, _ntp>, NumericalDiffMode::Central> *numDiff;
    LevenbergMarquardt<NumericalDiff<GnCost<T, _nx, _nu, _np, _ng, _ntp>, NumericalDiffMode::Central> > *lm;
/*#else 
    SampleNumericalDiff<GnCost<T, _nx, _nu, _np, _ng, _ntp> > *numDiff;
    LevenbergMarquardt<SampleNumericalDiff<GnCost<T, _nx, _nu, _np, _ng, _ntp> > > *lm;
#endif
    */
  };

  
// Generic functor
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  const int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

  // you should define that in the subclass :
//  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};


 template <typename T, 
   int _nx = Dynamic,
   int _nu = Dynamic,
   int _np = Dynamic,
   int _ng = Dynamic,
   int _ntp = Dynamic>
   struct GnCost : Functor<double> {  

   GnDocp<T, _nx, _nu, _np, _ng, _ntp> *docp;

   GnCost<T, _nx, _nu, _np, _ng, _ntp>(int inputs, int values) : Functor<double>(inputs, values), docp(0) {};
   
   typedef Matrix<double, _nx, 1> Vectornd;
   typedef Matrix<double, _nu, 1> Vectorcd;
   typedef Matrix<double, _ng, 1> Vectorgd;
   
   int operator()(const VectorXd &s, VectorXd &fvec) const
   {
     assert(docp);
     docp->tparam.From(docp->ts, docp->xs, docp->us, s, docp->p);

     //#DEBUG Number of function evaluations:
     //cout<<++(docp->nofevaluations)<<endl;
     ++(docp->nofevaluations);

     //cout<<"Docp->Xsback: "<<(docp->xs).back().transpose()<<endl;//what is this back doing?
     if(docp->p)
     {
       cout<<"docp->p: "<<(docp->p)->transpose()<<endl;
     }
     //cout<<"Xf: "<<docp->xs[(docp->xs).size()-1].transpose()<<endl;
     /*
       for (int i = 0, k = 0; k < docp.us.size(); ++k) {
       memcpy(docp.us[k].data(), s.data() + i, docp.sys.U.n*sizeof(double));
       i += docp.sys.U.n;
      }
      docp.Update(false);
      */
     
     Vectorgd &g = ((LsCost<T, _nx, _nu, _np, _ng>&)docp->cost).g;
     
     int i = 0;
      for (int k = 0; k < docp->us.size(); ++k) {
        const double &t = docp->ts[k];
        double h = docp->ts[k+1] - t;
        const T& x = docp->xs[k];
        const Vectorcd& u = docp->us[k];
        
        ((LqCost<T, _nx, _nu, _np, _ng>&)docp->cost).Res(g, t, x, u, h, docp->p);
        //double cost_debug = ((LqCost<T, _nx, _nu, _np, _ng>&)docp->cost).L(t, x, u, h, docp->p);
        memcpy(fvec.data() + i, g.data(), g.size()*sizeof(double));
        i += g.size();
        /*cout<<"docp->us["<<k<<"]: "<<(docp->us[k].transpose())<<endl;
        cout<<"g["<<k<<"]: "<<g.transpose()<<endl;
        cout<<"gcost["<<k<<"]: "<<0.5*(g.transpose()*g)<<endl;
        cout<<"Lcost["<<k<<"]: "<<cost_debug<<endl;
        // 
        ((LqCost<T, _nx, _nu>&)docp->cost).ResX(rx, t, x, h);
        memcpy(fvec.data() + i, rx.data(), rx.size()*sizeof(double));
        i += rx.size();

        ((LqCost<T, _nx, _nu>&)docp->cost).ResU(ru, t, u, h);
        memcpy(fvec.data() + i, ru.data(), ru.size()*sizeof(double));
        i += ru.size();
        */
      }
      
      ((LqCost<T, _nx, _nu, _np, _ng>&)docp->cost).Res(g, 
                                                       docp->ts.back(), 
                                                       docp->xs.back(), 
                                                       docp->us.back(), 0);
      //cout<<"g[end]: "<<g.transpose()<<endl;
      memcpy(fvec.data() + i, g.data(), g.size()*sizeof(double));
      docp->J = 0.5*(fvec.transpose()*fvec)(0);
      //cout<<"s: "<<s.transpose()<<endl;
      //cout<<"fvec: "<<fvec.transpose()<<endl;
      //cout<<"Cost: "<<(docp->J)<<endl;
      //getchar(); // #DEBUG
      return 0;
    }
  };



#ifdef GCOP_GNDOCP_CERES

  // A cost functor that implements the residual r = 10 - x.
  template <typename T, int _nx = Dynamic, int _nu = Dynamic>
    struct GnCost {  
      GnCost<T, _nx, _nu>(GnDocp<T, _nx, _nu> &docp) : docp(docp) {};

      GnDocp<T, _nx, _nu> &docp;

    typedef Matrix<double, _nx, 1> Vectornd;
    typedef Matrix<double, _nu, 1> Vectorcd;

    
    //   bool operator()(const double* const s, double* residual) const {
    bool operator()(double const* const* parameters, double* residuals) const {

      const double *s = parameters[0];
      for (int i = 0, k = 0; k < docp.us.size(); ++k) {
        memcpy(docp.us[k].data(), s + i, docp.sys.U.n*sizeof(double));
        i += docp.sys.U.n;
      }
      docp.Update(false);

      Vectorcd ru;
      Vectornd rx;
      int i = 0;
      for (int k = 0; k < docp.us.size(); ++k) {
        const double &t = docp.ts[k];
        double h = docp.ts[k+1] - t;
        const T& x = docp.xs[k];
        const Vectorcd& u = docp.us[k];

        ((LqCost<T, _nx, _nu>&)docp.cost).ResX(rx, t, x, h);
        memcpy(residuals + i, rx.data(), rx.size()*sizeof(double));
        i += rx.size();

        ((LqCost<T, _nx, _nu>&)docp.cost).ResU(ru, t, u, h);
        memcpy(residuals + i, ru.data(), ru.size()*sizeof(double));
        i += ru.size();
      }
      if (docp.ts.size() > 1) {
        double h = docp.ts.back() - docp.ts[docp.ts.size() - 2];
        ((LqCost<T, _nx, _nu>&)docp.cost).ResX(rx, docp.ts.back(), docp.xs.back(), h);
        memcpy(residuals + i, rx.data(), rx.size()*sizeof(double));
      }
      return true;
    }
  };
#endif

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int _nx, int _nu, int _np, int _ng, int _ntp> 
    GnDocp<T, _nx, _nu, _np, _ng, _ntp>::GnDocp(System<T, _nx, _nu, _np> &sys, 
                                                LqCost<T, _nx, _nu, _np, _ng> &cost,
                                                Tparam<T, _nx, _nu, _np, _ntp> &tparam,
                                                vector<double> &ts, 
                                                vector<T> &xs, 
                                                vector<Matrix<double, _nu, 1> > &us,
                                                Matrix<double, _np, 1> *p,
                                                bool update) : 
    Docp<T, _nx, _nu, _np>(sys, cost, ts, xs, us, p, false), tparam(tparam),
    inputs(tparam.ntp),
    values(cost.ng*xs.size()), s(inputs), 
    functor(0), numDiff(0), lm(0), 
//#ifndef USE_SAMPLE_NUMERICAL_DIFF
    numdiff_stepsize(1e-8)
//#else
    //numdiff_stepsize(1e-4)
//#endif
    {
      if(update)
        this->Update(false);//No need of derivatives

      ++(this->nofevaluations);

      tparam.To(s, this->ts, this->xs, this->us, this->p);

      cout <<"ntp=" <<tparam.ntp << endl;
    }
  
  template <typename T, int _nx, int _nu, int _np, int _ng, int _ntp> 
    GnDocp<T, _nx, _nu, _np, _ng, _ntp>::~GnDocp()
    {
      delete lm;
      delete numDiff;
      delete functor;
    }  

  template <typename T, int _nx, int _nu, int _np, int _ng, int _ntp> 
    void GnDocp<T, _nx, _nu, _np, _ng, _ntp>::Reset() {
      delete lm;
      lm = NULL;
      return;
    }
  
  template <typename T, int _nx, int _nu, int _np, int _ng, int _ntp> 
    void GnDocp<T, _nx, _nu, _np, _ng, _ntp>::Iterate() {

    if (!lm) {
      functor = new GnCost<T, _nx, _nu, _np, _ng, _ntp>(inputs, values);
      functor->docp = this;
//#ifndef USE_SAMPLE_NUMERICAL_DIFF
      numDiff = new NumericalDiff<GnCost<T, _nx, _nu, _np, _ng, _ntp>, NumericalDiffMode::Central>(*functor,numdiff_stepsize);
      lm = new LevenbergMarquardt<NumericalDiff<GnCost<T, _nx, _nu, _np, _ng, _ntp>, NumericalDiffMode::Central> >(*numDiff);
/*#else 
      numDiff = new SampleNumericalDiff<GnCost<T, _nx, _nu, _np, _ng, _ntp> >(*functor,numdiff_stepsize);
      lm = new LevenbergMarquardt<SampleNumericalDiff<GnCost<T, _nx, _nu, _np, _ng, _ntp> > >(*numDiff);
#endif
      */
      lm->parameters.maxfev = 1e6;//Maximum nof evaluations is very high
      cout<<"Initializing..."<<endl;
      info = lm->minimizeInit(s);
      cout <<"info="<<info <<endl;
    }

    /*
    for (int i = 0, k = 0; k < this->us.size(); ++k) {
      memcpy(s.data() + i, this->us[k].data(), this->sys.U.n*sizeof(double));
      i += this->sys.U.n;
    }
    */

    //    lm.parameters.maxfev=10000;
    info = lm->minimizeOneStep(s);
    
    cout <<"info="<<info <<endl;
    // check return values
    // VERIFY_IS_EQUAL(info, 1);
    //   VERIFY_IS_EQUAL(lm.nfev(), 26);
    
#ifdef GCOP_GNDOCP_CERES

    //  google::InitGoogleLogging(argv[0]);

  int npar = this->sys.U.n*this->us.size();
  int nres = ((LqCost<T, _nx, _nu, _np, _ng>&)cost).gn*this->xs.size();

  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  double s[npar];
  for (int i = 0, k = 0; k < this->us.size(); ++k) {
    memcpy(s + i, this->us[k].data(), this->sys.U.n*sizeof(double));
    i += this->sys.U.n;
  }

  //  const double initial_x = x;

  // Build the problem.
  Problem problem;

  //  CostFunction* cost_function =
  //    new NumericDiffCostFunction<CostFunctor, CENTRAL, nres, _npar> (new CostFunctor);
  //  problem.AddResidualBlock(cost_function, NULL, &x);
  
  DynamicNumericDiffCostFunction<GnCost<T,n,c> > *cost_function = new DynamicNumericDiffCostFunction<GnCost<T,n,c> >(new GnCost<T,n,c>(*this));
  cost_function->AddParameterBlock(npar);
  cost_function->SetNumResiduals(nres);

  problem.AddResidualBlock(cost_function, NULL, s);

  // Run the solver!
  Solver::Options options;

  options.linear_solver_type = ceres::DENSE_QR;
  //  options.num_threads = 8;
  options.minimizer_progress_to_stdout = true;
  options.max_num_iterations = 1;
  Solver::Summary summary;
  Solve(options, &problem, &summary);

  std::cout << summary.BriefReport() << "\n";
  //  std::cout << "x : " << initial_x
  #endif
    
  }
}

#endif
