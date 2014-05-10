#ifndef GCOP_GNDOCP_H
#define GCOP_GNDOCP_H

#include "docp.h"
#include "lqcost.h"

#include <unsupported/Eigen/NonLinearOptimization>

#ifdef GCOP_GNDOCP_CERES
#include "ceres/ceres.h"
#include "glog/logging.h"
#endif

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
 
  template <typename T, int n = Dynamic, int c = Dynamic> class GnDocp : 
    public Docp<T, n, c>{
    
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
    GnDocp(System<T, Vectorcd, n, c> &sys, LqCost<T, n, c> &cost,
           vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, bool update = true);
    
    virtual ~GnDocp();

    /**
     * Perform one DOCP iteration. Internally calls:
     * are updated. 
     */
    void Iterate();


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


  template <typename T, int n = Dynamic, int c = Dynamic>
    struct GnCost : Functor<double> {  
    GnCost<T,n,c>(GnDocp<T,n,c> &docp, int inputs, int values) : Functor<double>(inputs, values), docp(docp) {};
    
    GnDocp<T,n,c> &docp;
    
    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;

    int operator()(const VectorXd &s, VectorXd &fvec) const
    {
      for (int i = 0, k = 0; k < docp.us.size(); ++k) {
        memcpy(docp.us[k].data(), s.data() + i, docp.sys.c*sizeof(double));
        i += docp.sys.c;
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
        
        ((LqCost<T,n,c>&)docp.cost).ResX(rx, t, x, h);
        memcpy(fvec.data() + i, rx.data(), rx.size()*sizeof(double));
        i += rx.size();

        ((LqCost<T,n,c>&)docp.cost).ResU(ru, t, u, h);
        memcpy(fvec.data() + i, ru.data(), ru.size()*sizeof(double));
        i += ru.size();
      }
      if (docp.ts.size() > 1) {
        double h = docp.ts.back() - docp.ts[docp.ts.size() - 2];
        ((LqCost<T,n,c>&)docp.cost).ResX(rx, docp.ts.back(), docp.xs.back(), h);
        memcpy(fvec.data() + i, rx.data(), rx.size()*sizeof(double));
      }
      return 0;
    }
  };


#ifdef GCOP_GNDOCP_CERES

  // A cost functor that implements the residual r = 10 - x.
  template <typename T, int n = Dynamic, int c = Dynamic>
    struct GnCost {  
      GnCost<T,n,c>(GnDocp<T,n,c> &docp) : docp(docp) {};

      GnDocp<T,n,c> &docp;

    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;

    
    //   bool operator()(const double* const s, double* residual) const {
    bool operator()(double const* const* parameters, double* residuals) const {

      const double *s = parameters[0];
      for (int i = 0, k = 0; k < docp.us.size(); ++k) {
        memcpy(docp.us[k].data(), s + i, docp.sys.c*sizeof(double));
        i += docp.sys.c;
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

        ((LqCost<T,n,c>&)docp.cost).ResX(rx, t, x, h);
        memcpy(residuals + i, rx.data(), rx.size()*sizeof(double));
        i += rx.size();

        ((LqCost<T,n,c>&)docp.cost).ResU(ru, t, u, h);
        memcpy(residuals + i, ru.data(), ru.size()*sizeof(double));
        i += ru.size();
      }
      if (docp.ts.size() > 1) {
        double h = docp.ts.back() - docp.ts[docp.ts.size() - 2];
        ((LqCost<T,n,c>&)docp.cost).ResX(rx, docp.ts.back(), docp.xs.back(), h);
        memcpy(residuals + i, rx.data(), rx.size()*sizeof(double));
      }
      return true;
    }
  };
#endif

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int n, int c> 
    GnDocp<T, n, c>::GnDocp(System<T, Matrix<double, c, 1>, n, c> &sys, 
                            LqCost<T, n, c> &cost,
                            vector<double> &ts, 
                            vector<T> &xs, 
                            vector<Matrix<double, c, 1> > &us,
                            bool update) : 
    Docp<T,n,c>(sys, cost, ts, xs, us, update) {
  }
  
  template <typename T, int n, int c> 
    GnDocp<T, n, c>::~GnDocp()
    {
    }  
  
  template <typename T, int n, int c> 
    void GnDocp<T, n, c>::Iterate() {
    
    int info;
    double fnorm, covfac;

    int inputs = this->sys.c*this->us.size();
    int values = this->sys.c*this->us.size() + this->sys.n*this->xs.size();
    
    VectorXd s(inputs);

    for (int i = 0, k = 0; k < this->us.size(); ++k) {
      memcpy(s.data() + i, this->us[k].data(), this->sys.c*sizeof(double));
      i += this->sys.c;
    }    

    // do the computation
    GnCost<T,n,c> functor(*this, inputs, values);
    
    NumericalDiff<GnCost<T,n,c> > numDiff(functor);

    LevenbergMarquardt<NumericalDiff<GnCost<T,n,c> > > lm(numDiff);

    //    lm.parameters.maxfev=10000;
    info = lm.minimize(s);

    cout <<"info="<<info <<endl;
    // check return values
    // VERIFY_IS_EQUAL(info, 1);
    //   VERIFY_IS_EQUAL(lm.nfev(), 26);
    
#ifdef GCOP_GNDOCP_CERES

    //  google::InitGoogleLogging(argv[0]);

  int npar = this->sys.c*this->us.size();
  int nres = this->sys.c*this->us.size() + this->sys.n*this->xs.size();

  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  double s[npar];
  for (int i = 0, k = 0; k < this->us.size(); ++k) {
    memcpy(s + i, this->us[k].data(), this->sys.c*sizeof(double));
    i += this->sys.c;
  }

  //  const double initial_x = x;

  // Build the problem.
  Problem problem;

  //  CostFunction* cost_function =
  //    new NumericDiffCostFunction<CostFunctor, CENTRAL, nres, npar> (new CostFunctor);
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
