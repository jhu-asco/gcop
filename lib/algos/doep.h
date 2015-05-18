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
     * @param ts1 (N1+1) sequence of discrete times for sensor measurements
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    Doep(System<T, nx, nu, np> &sys, Sensor<T1, nx1, nu, np, Tz, nz> &sensor, 
         SensorCost<T, nx, nu, np, Tz, nz> &cost, 
         vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
         Vectormd &p, vector<double> &ts1, Func_type _project=NULL, bool update = true);
    
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

    std::vector<double> &ts1; ///< times (N+1) vector for sensor measurements

    std::vector<T> &xs;      ///< states (N+1) vector

    std::vector<Vectorcd> &us;      ///< controls (N) vector
    
    Vectormd &p;               ///< parameter vector

    std::vector<Vectornd> ws; ///< Internal process noise

    std::vector<Vectorrd> zs; ///< Internal sensor measurements
    
    bool debug;  ///< whether to display debugging info
    
    double eps;  ///< epsilon used for finite differences

    double J;    ///< Optimal cost at any point of time

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
                              vector<double> &ts1, Func_type _project, bool update) : 
    sys(sys), sensor(sensor), cost(cost), ts(ts), xs(xs), us(us), p(p), ts1(ts1), project(_project),
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
      zs.resize(ts1.size()); //Internal sensor measurements
      if(nz == Dynamic)
      {
        for(int i = 0; i<ts1.size(); ++i)
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

    sys.Reset(xs[0],ts[0]);//Reset
    int sensor_index = 0;
    assert(ts1[0]>= ts[0]);//Sensor measurements are after the control started
    assert(ts1.back() <= ts.back());//Sensor measurements are before control ended

    for (int k = 0; k < N; ++k) {
      double h = ts[k+1] - ts[k];
      sys.Step(xs[k+1], us[k], h, ws[k], &p);
      //cout<<"ts1["<<sensor_index<<"]: "<<ts1[sensor_index]<<"\t"<<ts[k]<<"\t"<<ts[k+1]<<endl;//#DEBUG
      while((ts1[sensor_index] - ts[k])>= 0 && (ts1[sensor_index] - ts[k+1]) < 0)//Nearest state to find the sensor measurement
      {
        //#TODO Add Linear Interpolation if necessary 
        int near_index = (ts1[sensor_index] - ts[k]) > -(ts1[sensor_index] - ts[k+1])?(k+1):k;
        //Project the state
        if(std::is_same<T, T1>::value)
        {
          x1 = (T1)(*((T1 *)(&xs[near_index])));//Hacky way of copying pointer over
          sensor(zs[sensor_index], ts[near_index], x1, us[near_index], &p);
        }
        else
        {
          assert(project != NULL);// Can be removed
          project(xs[near_index],x1);
          sensor(zs[sensor_index], ts[near_index], x1, us[near_index], &p);
          //cout<<"UPDATE: "<<"Zs["<<sensor_index<<"]: "<<zs[sensor_index].transpose()<<endl;//#DEBUG
        }
        if(sensor_index == (ts1.size()-1))
          break;
        sensor_index = sensor_index < (ts1.size()-1)?sensor_index+1:sensor_index;
      }
      //getchar();//#DEBUG
    }
    assert(sensor_index == (ts1.size()-1));//Assert that we collected all the sensor data
  }  

  template <typename T, int nx, int nu, int np, typename Tz, int nz, typename T1, int nx1> 
    void Doep<T, nx, nu, np, Tz, nz, T1, nx1>::Iterate() {
    cout << "[W] Doep::Iterate: subclasses should implement this!" << endl;
  }
}

#endif
