#ifndef GCOP_DOEP_H
#define GCOP_DOEP_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include "system.h"
#include "sensor.h"
#include "sensorcost.h"
#include <cmath>
#include "rn.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic,
    typename Tz=VectorXd,
    int nz = Dynamic,
    typename T1=T,
    int nx1=nx> class Doep {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, np, 1> Vectormd;
    typedef Matrix<double, nz, 1> Vectorrd;

    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  
    
    typedef void(*Func_type)(const T &, T1 &);

  public:
    /**
     * Create a discrete optimal estimation problem using a system, sensor, a cost, and 
     * a trajectory given by a sequence of times, states, and controls. 
     * The controls us are provided. The times ts must be given, the initial state xs[0] must be set.
     *
     * After initialization, every call to Iterate() will optimize the parameters p and  
     * process noise wi are varied to fit the measurements in the sensorcost
     *
     * 
     * @param sys system
     * @param sensor sensor
     * @param cost cost
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    Doep(System<T, nx, nu, np> &sys, Sensor<T1, nx1, nu, np, Tz, nz> &sensor, 
         SensorCost<T, nx, nu, np, Tz, nz> &cost, 
         vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
         Vectormd &p, Func_type _project=NULL, bool update = true);
    
    virtual ~Doep();

    /**
     * Perform one DOCP iteration. 
     */
    virtual void Iterate();

    /**
     * Update the trajectory and  sensor measurements 
     * @param der NOT USED RIGHT NOW
     */
    void Update(bool der = true);
    
    System<T, nx, nu, np> &sys;    ///< dynamical system 

    Sensor<T1, nx1, nu, np, Tz, nz> &sensor;    ///< sensor

    SensorCost<T, nx, nu, np, Tz, nz> &cost;     ///< given cost function

    std::vector<double> &ts; ///< times (N+1) vector

    std::vector<T> &xs;      ///< states (N+1) vector

    std::vector<Vectorcd> &us;      ///< controls (N) vector
    
    Vectormd &p;               ///< parameter vector

    std::vector<Vectornd> ws; ///< Internal process noise

    std::vector<Vectorrd> zs; ///< Internal sensor measurements
    
    bool debug;  ///< whether to display debugging info
    
    double eps;  ///< epsilon used for finite differences

    Func_type project; ///<Project T to T1 will add jacobians later

  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int nx, int nu, int np, typename Tz, int nz, typename T1, int nx1> 
    Doep<T, nx, nu, np, Tz, nz, T1, nx1>::Doep(System<T, nx, nu, np> &sys, Sensor<T1, nx1, nu, np, Tz, nz> &sensor,
                              SensorCost<T, nx, nu, np, Tz, nz> &cost, 
                              vector<double> &ts, 
                              vector<T> &xs, 
                              vector<Matrix<double, nu, 1> > &us,
                              Matrix<double, np, 1> &p,
                              Func_type _project, bool update) : 
    sys(sys), sensor(sensor), cost(cost), ts(ts), xs(xs), us(us), p(p), project(_project),
    debug(true), eps(1e-3)
    {
      int N = us.size();
      assert(N > 0);
      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      ws.resize(N);//Internal process noise
      if(nx == Dynamic)
      {
        for(int i =0; i<N; ++i)
        {
          ws[i].resize(sys.X.n);
        }
      }
      zs.resize(N); //Internal sensor measurements
      if(nz == Dynamic)
      {
        for(int i = 0; i<N; ++i)
        {
          zs[i].resize(sensor.Z.n);
        }
      }
      for(int i =0; i<N; ++i)
      {
        ws[i].setZero();//Set the noise to zero when constructing
      }

      if (update) {
        Update();
      }
    }
  
  template <typename T, int nx, int nu, int np, typename Tz, int nz, typename T1, int nx1> 
    Doep<T, nx, nu, np, Tz, nz, T1, nx1>::~Doep()
    {
    }
  
  template <typename T, int nx, int nu, int np, typename Tz, int nz, typename T1, int nx1> 
    void Doep<T, nx, nu, np, Tz, nz, T1, nx1>::Update(bool der) {

    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nz, 1> Vectorrd;
    
   
    int N = us.size();
    T1 x1; //Temporary projection for finding sensor measurements

    sys.reset(xs[0],ts[0]);//Reset

    for (int k = 0; k < N; ++k) {
      double h = ts[k+1] - ts[k];
      sys.Step3(xs[k+1], us[k], ws[k], h, &p);
      //Project the state
      if(std::is_same<T, T1>::value)
      {
        x1 = (T1)(*((T1 *)(&xs[k])));//Hacky way of copying pointer over
        sensor(zs[k], ts[k], x1, us[k], &p);
      }
      else
      {
        assert(project != NULL);// Can be removed
        project(xs[k],x1);
        sensor(zs[k], ts[k], x1, us[k], &p);
        //cout<<"Zs["<<k<<"]: "<<zs[k].transpose()<<endl;//#DEBUG
      }
    }
  }  

  template <typename T, int nx, int nu, int np, typename Tz, int nz, typename T1, int nx1> 
    void Doep<T, nx, nu, np, Tz, nz, T1, nx1>::Iterate() {
    cout << "[W] Doep::Iterate: subclasses should implement this!" << endl;
  }
}

#endif
