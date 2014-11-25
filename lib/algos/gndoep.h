#ifndef GCOP_GNDOEP_H
#define GCOP_GNDOEP_H

#include "doep.h"
#include "lqsensorcost.h"

#include <unsupported/Eigen/NonLinearOptimization>

namespace gcop {

 
  using namespace std;
  using namespace Eigen;

 
  template <typename T, 
    int _nx, 
    int _nu,
    int _np,
    int _ng,
    typename Tz,
    int _nz,
    typename T1,
    int _nx1> class GnCost;
  
  template <typename T = VectorXd, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic,
    int _ng = Dynamic,
    typename Tz = VectorXd,
    int _nz = Dynamic,
    typename T1 = T,
    int _nx1 = _nx> 
    class GnDoep : public Doep<T, _nx, _nu, _np, Tz, _nz, T1, _nx1> {
    
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

    typedef Matrix<double, _nz, 1> Vectorrd;
    typedef Matrix<double, _nz, _nz> Matrixrd;
    typedef Matrix<double, _nz, _nx> Matrixrnd;
    typedef Matrix<double, _nz, _nu> Matrixrcd;
    typedef Matrix<double, _nz, _np> Matrixrmd;

    typedef void(*Func_type)(const T &, T1 &);

  public:
    /**
     * Create an optimal estimation problem using a system, sensor, a cost, and
     * sensor measurements and parameters of the system.
     * 
     * The optimal estimation is done over selecting the right parameters for the system
     * and also process noise to ensure the estimation cost on sensor measurements and noise is minimized
     * Parameters can also have a prior which is added into the cost function. The parameters passed in here
     * are used as initial guess to the system
     *
     * You should set the ts[1:N+1], xs[0]; Known controls us
     *
     * After initialization, every call to Iterate() will optimize the 
     * parameters of the system, process noise and modify them accordingly.
     *
     * @param sys system
     * @param cost cost
     * @param tparam trajectory parametrization
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param p static parameter vector initial guess of it
     * @param ng 
     * @param _project Optional function to project T into T1. Not needed if T == T1
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */


    GnDoep(System<T, _nx, _nu, _np> &sys, Sensor<T1, _nx1, _nu, _np, Tz, _nz> &sensor,
           LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz> &cost, 
           vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, Vectormd &p,
           Func_type _project=NULL, bool update = true);
    
    virtual ~GnDoep();

    /**
     * Perform one DOCP iteration. Internally calls:
     * are updated. 
     */
    void Iterate();

    int info;
    double fnorm, covfac;

    int inputs;
    int values;
    
    VectorXd s;  ///< optimization vector 
    
    GnCost<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1> *functor;
    NumericalDiff<GnCost<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1> > *numDiff;
    LevenbergMarquardt<NumericalDiff<GnCost<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1> > > *lm;
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
   typename Tz=VectorXd,
   int _nz = Dynamic,
   typename T1 = T,
   int _nx1 = _nx>
   struct GnCost : Functor<double> {  
   GnCost<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1>(int inputs, int values) : Functor<double>(inputs, values) {};
   
   GnDoep<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1> *doep;
   
   typedef Matrix<double, _nx, 1> Vectornd;
   typedef Matrix<double, _np, 1> Vectormd;
   typedef Matrix<double, _nu, 1> Vectorcd;
   typedef Matrix<double, _ng, 1> Vectorgd;
   typedef Matrix<double, _nz, 1> Vectorrd;
   
   int operator()(const VectorXd &s, VectorXd &fvec) const
   {
     assert(doep);

     const int &np = doep->sys.P.n;
     const int &nz = doep->sensor.Z.n;
     const int &nw = doep->sys.X.n;
     const int N = doep->us.size();

     //Update inputs: process noise ws[i], p:
     doep->p = s.tail(np);
     for(int i = 0; i < N; i++)
     {
       doep->ws[i] = s.segment(i*nw,nw);
     }
     //Update sensor measurements using the input parameters and noise ws:
     doep->Update(false);

     Vectorgd &g = ((LsSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>&)doep->cost).g;
     Vectormd &gp = ((LsSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>&)doep->cost).gp;
     fvec.tail(np).setZero();

     //Compute and set residuals
     int i = 0;
     for(int k = 0; k< N; ++k)
     {
       double h = doep->ts[k+1] - doep->ts[k];
       ((LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>&)doep->cost).Res(g, doep->ts[k], doep->zs[k], doep->ws[k], doep->p, h);
       memcpy(fvec.data() + i, g.data(), (nw+nz)*sizeof(double));
       //cout<<"Res: "<<g<<endl;
       //cout<<"i: "<<i<<endl;
       i += (nw+nz);
       fvec.tail(np) += g.tail(np);//The tail is a constant residual for parameters
     }
     ((LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>&)doep->cost).Resp(gp, doep->p);
     fvec.tail(np) += gp;
     //cout<<"Gp: "<<gp<<endl;

     if(doep->debug)
     {
       std::cout<<"Input: "<<s.transpose()<<endl;
       std::cout<<"Fvec: "<<fvec.transpose()<<endl;
       //std::cout<<"Resp: "<<fvec.tail(np).transpose()<<"\t"<<doep->p<<endl;
     }
       std::cout<<"Cost: "<<0.5*(fvec.transpose()*fvec)<<endl;
     return 0;
   }
  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz, typename T1, int _nx1> 
    GnDoep<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1>::GnDoep(System<T, _nx, _nu, _np> &sys, Sensor<T1, _nx1, _nu, _np, Tz, _nz> &sensor,
                                                LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz> &cost,
                                                vector<double> &ts, 
                                                vector<T> &xs, 
                                                vector<Vectorcd > &us,
                                                Vectormd &p,
                                                Func_type _project,
                                                bool update) : 
    Doep<T, _nx, _nu, _np, Tz, _nz, T1, _nx1>(sys, sensor, cost, ts, xs, us, p, _project, update),
    inputs(us.size()*sys.X.n + sys.P.n),
    values((sys.X.n + sensor.Z.n)*us.size()+sys.P.n), s(inputs), 
    functor(0), numDiff(0), lm(0)
    {
      cout <<"inputs=" <<inputs<<" values= "<<values<< endl;
    }
  
  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz, typename T1, int _nx1> 
    GnDoep<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1>::~GnDoep()
    {
      delete lm;
      delete numDiff;
      delete functor;
    }  
  
  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz, typename T1, int _nx1> 
    void GnDoep<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1>::Iterate() {

    if (!lm) {
      functor = new GnCost<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1>(inputs, values);
      functor->doep = this;
      numDiff = new NumericalDiff<GnCost<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1> >(*functor);
      lm = new LevenbergMarquardt<NumericalDiff<GnCost<T, _nx, _nu, _np, _ng, Tz, _nz, T1, _nx1> > >(*numDiff);

      s.setZero();//Set the initial noise ws to zero
      const int &np = this->sys.P.n;
      s.tail(np) = this->p;//Set the system parameters to initial guess
    }

    /*
    for (int i = 0, k = 0; k < this->us.size(); ++k) {
      memcpy(s.data() + i, this->us[k].data(), this->sys.U.n*sizeof(double));
      i += this->sys.U.n;
    }
    */

    //    lm.parameters.maxfev=10000;
    info = lm->minimize(s);
    
    cout <<"info="<<info <<endl;
    // check return values
    // VERIFY_IS_EQUAL(info, 1);
    //   VERIFY_IS_EQUAL(lm.nfev(), 26);
  }
}

#endif
