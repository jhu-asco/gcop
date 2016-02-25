#ifndef GCOP_QROTORID_GNDOCP_H
#define GCOP_QROTORID_GNDOCP_H

#include "docp.h"
#include "lqcost.h"
#include "tparam.h"
#include "qrotoridmodel.h"
#include <fstream>

#include <unsupported/Eigen/NonLinearOptimization>

namespace gcop {

  using namespace std;
  using namespace Eigen;

  struct GnCost;//Forward Declaration

  struct Obstacle{
      double radius;//Cylinder
      Vector3d center;//< Center of cylinder
      Vector3d axis;///< Axis of cylinder (If axis is zero, it becomes a sphere)
      int id;///< Cylinder: 0, Sphere 1, Plane 2
      Obstacle():radius(0.1), center(0,0,0), axis(0,0,1), id(0){

      }
      void set(VectorXd info, int id_ = 0)
      {
        id = id_;
        if(id == 2)
          radius = 0;
        else
          radius = info(0);
        center = info.segment<3>(1);
        axis = info.segment<3>(4);
      }
      /**
       * @brief Distance  between the current position and Obstacle
       * @param pos Current Quad Position
       * @param invcov
       * @return
       */
      void Distance(const Vector3d &pos, const Matrix3d &invcov, double &ko, double *res_ref)
      {
          Vector3d error = center - pos;
          if(id == 0 || id == 1)
            error = error - (error.dot(axis))*axis;//Remove component along axis
          else
            error = (error.dot(axis))*axis;//Only care about component along axis
          double scale_ellipsoid;

          scale_ellipsoid = 0.5*sqrt(error.transpose()*invcov*error);//2 Stdeviations
          if(std::isinf(scale_ellipsoid) || !isfinite(scale_ellipsoid))
          {
            cout<<"Received Nan or Inf on scale ellipsoid"<<endl;
            scale_ellipsoid = 1e10;
          }
          //scale_ellipsoid= sqrt((error.cwiseQuotient(2*stdev).array().square()).sum());
          //scale_ellipsoid = error.norm()/(2*stdev.maxCoeff());//Consider this as a sphere conservative estimate

          double scale_radius = radius/error.norm();
          double scale_error = 1 - 1.0/scale_ellipsoid - scale_radius;
          Map<Vector3d> dist_vec(res_ref);
          if(scale_error > -1e-10)//Not penetrating obstacle
          {
            dist_vec.setZero();
          }
          else
          {
            //cout<<"Colliding: "<<pos.transpose()<<"\t"<<stdev.transpose()<<"\t"<<center.transpose()<<endl;
            dist_vec = ko*scale_error*error;
          }
      }
  };
 
    class QRotorIdGnDocp :
    public Docp<QRotorIDState, 15, 4, 13> {
    
    typedef Matrix<double, 19, 1> Vectorgd;
    typedef Matrix<double, 19, 15> Matrixgxd;
    typedef Matrix<double, 19, 4> Matrixgud;
    typedef Matrix<double, 19, 13> Matrixgpd;
    
    typedef Matrix<double, 15, 1> Vectornd;
    typedef Matrix<double, 4, 1> Vectorcd;

    typedef Matrix<double, 15, 15> Matrixnd;
    typedef Matrix<double, 4, 4> Matrixncd;
    typedef Matrix<double, 4, 15> Matrixcnd;
    typedef Matrix<double, 4, 15> Matrixcd;
    
    typedef Matrix<double, 13, 1> Vectormd;
    
    typedef Matrix<double, 13, 13> Matrixmd;
    typedef Matrix<double, 15, 13> Matrixnmd;
    typedef Matrix<double, 13, 15> Matrixmnd;
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
    QRotorIdGnDocp(System<QRotorIDState, 15, 4, 13> &sys,
           LsCost<QRotorIDState, 15, 4, 13, 19> &cost, Tparam<QRotorIDState, 15, 4, 13> &tparam,
           vector<double> &ts, vector<QRotorIDState> &xs, vector<Vectorcd> &us,
           Vectormd *p = 0, bool update = true);
    
    virtual ~QRotorIdGnDocp();

    /**
     * Perform one DOCP iteration. Internally calls:
     * are updated.
     */
    void Iterate();

    /**
     * @brief Add an obstacle to the gndocp problem
     * @param obs
     */
    void AddObstacles(vector<Obstacle> &obs);

    /**
     * @brief GenerateStdev : Generates xs mean and stdeviation from sampling the uncertainity in initial state and parameters
     * @param alpha Scaling parameter on how far the samples from mean
     * @param kappa Secondary scaling parameter usually set to zero
     * @param beta  Prior Information on distribution. Optimal value for gaussian is 2
     */
    void GenerateStdev(double alpha = 1e-3, double kappa = 0, double beta = 2);

    void Reset();
    Tparam<QRotorIDState, 15, 4, 13> &tparam;

    int info;
    double fnorm, covfac;

    int inputs;
    int values;

    double numdiff_stepsize;   ///< The step size for perturbations
     
    VectorXd s;  ///< optimization vector
    
    GnCost *functor;

    Matrixmd stdev_params;///< Covariance of parameters for sampling

    Matrixnd stdev_initial_state;///< Covariance of initial state

    vector<Matrix3d> xs_invcov;///< Inverse of Covariance of posn from sampling for Obstacle detection

    Matrix<double,3,Dynamic> sample_mean_params;///< Sample of mean trajectory

    double ko;///< Obstacle Avoidance Gain

    int number_obstacles;///< Number of obstacles

    protected:

    vector<Obstacle> obstacles;///< vector of obstacles

    MatrixXd sample_xs;///< Temporary variable for holding sampled trajectories in Unscented Transform function

    NumericalDiff<GnCost, NumericalDiffMode::Central> *numDiff;

    LevenbergMarquardt<NumericalDiff<GnCost, NumericalDiffMode::Central> > *lm;

    friend struct  GnCost;
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

};


   struct GnCost : Functor<double> {

   QRotorIdGnDocp *docp;

   GnCost(int inputs, int values) : Functor<double>(inputs, values), docp(0) {}
   
   typedef Matrix<double, 15, 1> Vectornd;
   typedef Matrix<double, 4, 1> Vectorcd;
   typedef Matrix<double, 19, 1> Vectorgd;
   
   int operator()(const VectorXd &s, VectorXd &fvec) const
   {
     assert(docp);
     docp->tparam.From(docp->ts, docp->xs, docp->us, s, docp->p);

     //Evaluate Mean and Covariance of trajectories by sampling from dist of params and inputs us

     ++(docp->nofevaluations);

     Vectorgd &g = ((LsCost<QRotorIDState, 15, 4, 13, 19>&)docp->cost).g;
     

     //if(number_obstacles > 0)
       //docp->GenerateStdev();//Generate Mean and Stdev of trajectory

     double *fvec_ref = fvec.data();
     //Can parallelize this
      for (int k = 0; k < docp->us.size(); ++k) {
        const double &t = docp->ts[k];
        double h = docp->ts[k+1] - t;
        const QRotorIDState& x = docp->xs[k];
        const Vectorcd& u = docp->us[k];
        
         //Evaluate Regular Cost:
        ((LsCost<QRotorIDState, 15, 4, 13, 19>&)docp->cost).Res(g, t, x, u, h, docp->p);
        memcpy(fvec_ref, g.data(), g.size()*sizeof(double));
        fvec_ref = fvec_ref + g.size();//Forward reference
        //Evaluate Cost from Obstacles:
        for(int obs_ind = 0;obs_ind < docp->number_obstacles; obs_ind++)
        {
            //const Vector3d &pos = docp->sample_mean_params.col(k);
            const Vector3d &pos = x.p;
            const Matrix3d &std = docp->xs_invcov.at(k);
            docp->obstacles[obs_ind].Distance(pos, std, docp->ko, fvec_ref);
            fvec_ref = fvec_ref+3;
        }
      }
      
      ((LsCost<QRotorIDState, 15, 4, 13, 19>&)docp->cost).Res(g,
                                                       docp->ts.back(),
                                                       docp->xs.back(), 
                                                       docp->us.back(), 0, docp->p);
      memcpy(fvec_ref, g.data(), g.size()*sizeof(double));
      fvec_ref = fvec_ref + g.size();//Forward reference

      //Evaluate Cost from Obstacles for final state
      for(int obs_ind = 0;obs_ind < docp->number_obstacles; obs_ind++)
      {
        //const Vector3d &pos = docp->sample_mean_params.col(docp->us.size());
        const Vector3d &pos = docp->xs.back().p;
        const Matrix3d &std = docp->xs_invcov.back();
        docp->obstacles[obs_ind].Distance(pos, std, docp->ko, fvec_ref);
        fvec_ref = fvec_ref + 3;
      }

      docp->J = 0.5*(fvec.transpose()*fvec)(0);
      //cout<<"s: "<<s.transpose()<<endl;
      //cout<<"fvec: "<<fvec.transpose()<<endl;
      //cout<<"Cost: "<<(docp->J)<<endl;
      //getchar(); // #DEBUG
      return 0;
    }
  };



  using namespace std;
  using namespace Eigen;
  
    QRotorIdGnDocp::QRotorIdGnDocp(System<QRotorIDState, 15, 4, 13> &sys,
                                                LsCost<QRotorIDState, 15, 4, 13, 19> &cost,
                                                Tparam<QRotorIDState, 15, 4, 13> &tparam,
                                                vector<double> &ts,
                                                vector<QRotorIDState> &xs,
                                                vector<Matrix<double, 4, 1> > &us,
                                                Matrix<double, 13, 1> *p,
                                                bool update) :
    Docp<QRotorIDState, 15, 4, 13>(sys, cost, ts, xs, us, p, false), tparam(tparam),
    inputs(tparam.ntp),
    values((cost.ng)*xs.size()), s(inputs),//Values is increased by one more to add obstacle avoidance term
    functor(0), numDiff(0), lm(0),
    numdiff_stepsize(1e-8), ko(0.01)
    {
      ++(this->nofevaluations);

      tparam.To(s, this->ts, this->xs, this->us, this->p);

      if(update)
        this->Update(false);//No need of derivatives

      cout <<"ntp=" <<tparam.ntp << endl;

      sample_mean_params.resize(3,xs.size());
      xs_invcov.resize(xs.size());
      stdev_initial_state.setZero();
      stdev_params.setZero();
    }
  
    QRotorIdGnDocp::~QRotorIdGnDocp()
    {
      delete lm;
      delete numDiff;
      delete functor;
    }  

    void QRotorIdGnDocp::Reset() {
      //if(ko >= 64)//Only reset if hitting bound
        //ko = 1.0;
      delete lm;
      lm = NULL;
      return;
    }

    void QRotorIdGnDocp::GenerateStdev(double alpha, double kappa, double beta)
    {
        QRotorIDState mean_initial_state = xs[0];
        Vectormd mean_params = *(this->p);
        int L = 28;//size of all parameters combined 15+13
        double lambda = (alpha*alpha)*(L+kappa)-L;
        double scales[4] = {sqrt(L+lambda), lambda/(L+lambda), (1/(2*(L+lambda))), (lambda/(L+lambda)) + (1 - alpha*alpha + beta)};
        int N = xs.size();//Number of states
        sample_xs.resize(3*N, 2*L);//Stdev from Unscented transform
        //Generate 2*L samples
        //Perturb in two directions
        int index_traj = 0;
        for(int k = -1; k < 2; k+=2)
        {
          for(int i = 0; i < 15;i++)
          {
            sys.X.Retract(xs[0],mean_initial_state,k*scales[0]*stdev_initial_state.col(i));//Change initial state
            *(this->p) = mean_params;
            this->Update(false);//Update states using current control
            for(int j = 0; j < N; j++)
            {
              sample_xs.block<3,1>(3*j,index_traj) = xs[j].p;
            }
            index_traj++;
          }
          for(int i = 0; i < 13;i++)
          {
            xs[0] = mean_initial_state;
            *(this->p) = mean_params + k*scales[0]*stdev_params.col(i);//Change input params
            this->Update(false);//Update states using current control and parameters
            for(int j = 0; j < N; j++)
            {
              sample_xs.block<3,1>(3*j,index_traj) = xs[j].p;
            }
            index_traj++;
          }
        }
        //Generate Trajectory without perturbations:
        xs[0] = mean_initial_state;
        *(this->p) = mean_params;
        this->Update(false);//Update states using current control and parameters
#pragma omp parallel for
        for(int j = 0; j < N; j++)
        {
          sample_mean_params.col(j) = xs[j].p;
        }
        
        Map<VectorXd> sample_mean_params_map(sample_mean_params.data(),3*N);
        //Trajectory Mean:
        VectorXd sample_mean = scales[2]*sample_xs.rowwise().sum() + scales[1]*sample_mean_params_map;
        //Center Samples:
        sample_xs = sample_xs.colwise() - sample_mean;
#pragma omp parallel for
        for(int j = 0; j < N; j++)
        {
          xs_invcov[j] = (scales[2]*sample_xs.middleRows<3>(3*j)*sample_xs.middleRows<3>(3*j).transpose()+ scales[3]*(sample_mean_params.col(j)-sample_mean.segment<3>(3*j))*(sample_mean_params.col(j)-sample_mean.segment<3>(3*j)).transpose()).selfadjointView<Upper>() ;
          xs_invcov[j] = xs_invcov[j].inverse().eval();//Invert covariance
        }
        //Map<VectorXd> xs_std_map(xs_std.data(),3*N);//xs_std_map shares data with xs_std
        //xs_std_map = (scales[2]*(((sample_xs.colwise() - sample_mean).array().square()).rowwise().sum())  + scales[3]*((sample_mean_params_map - sample_mean).array().square()) ).array().sqrt();//Stdeviation of xs
        //DEBUG:
        //cout<<"Xs_std: "<<xs_std.col(0).transpose()<<endl;
    }

    void QRotorIdGnDocp::AddObstacles(vector<Obstacle> &obs)
    {
        values = values +3*(obs.size() - obstacles.size())*(this->xs.size());
        obstacles = obs;//Copy Obstacles
        number_obstacles = obstacles.size();
    }
  
    void QRotorIdGnDocp::Iterate() {

    //cout<<"Obs Gain: "<<ko<<endl;//DEBUG

    if(obstacles.size() > 0)
      GenerateStdev();//Generate Mean and Stdev of trajectory

    if (!lm) {
      functor = new GnCost(inputs, values);
      functor->docp = this;
      numDiff = new NumericalDiff<GnCost, NumericalDiffMode::Central>(*functor,numdiff_stepsize);
      lm = new LevenbergMarquardt<NumericalDiff<GnCost, NumericalDiffMode::Central> >(*numDiff);

      //lm->parameters.maxfev = 5000;//Maximum nof evaluations is very high
      lm->parameters.maxfev = 2000;//Maximum nof evaluations is very high
      //cout<<"Initializing..."<<endl;
      info = lm->minimizeInit(s);
    }

    //tparam.From(ts,xs,us,s,p);
    for(int count = 0; count < 200; count++)
    {
      info = lm->minimizeOneStep(s);
      if(info != -1 && info != -2)
        break;
      if(obstacles.size() > 0)
        GenerateStdev();//Generate Mean and Stdev of trajectory
    }

    cout <<"info="<<info <<endl;
  }
}

#endif
