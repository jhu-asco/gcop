#ifndef GCOP_TRAJECTORY_H
#define GCOP_TRAJECTORY_H

#include <iostream>
#include <cstring>
#include <map>
#include "system.h"

namespace gcop {

  /**
   * Basic control system trajectory represented by a set of discrete nodes
   * 
   * This should serve as a base class for implementing trajectories for a particular
   * mechanical system described by the System class
   *
   * Author: Marin Kobilarov -- Copyright (C) 2006
   */
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic,
    int _ns = Dynamic> 
    class Trajectory {
  public:

    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  

    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, _ns, 1> Vectorsd;


  /**
   * Initialize a trajectory for the system sys
   * Once a trajectory is initialized with this constructor
   * one should call Init(sn,...) with the appropriate number of
   * discrete points desired. Otherwise the trajectory would contain
   * only one point (usually the origin)
   * @param sys mechanical system
   * @param sn trajectory parametrization vector length
   * @param N number of discrete segments (optional)
   */
  Trajectory(const System<T, nx, nu, np>& sys, int ns = 0, int N = 0);
  
  /**
   * Copy constructor
   * @param traj a trajectory
   */
  //  Trajectory(const Trajectory<T, nx, nu, np, _ns> &traj);
  
  /**
   * Copy assignment
   * @param traj a trajectory 
   */
  //  Trajectory& operator=(const Trajectory<T, nx, nu, np, _ns> &traj);
  
  
  /**
   * Initialize a trajectory from a serialized stream
   * @param sys mechanical system
   * @param istr stream
   */    
  //  Trajectory(const System<T, nx, nu, np>& sys, std::istream& istr);
  
  
  virtual ~Trajectory();
  
    /**
     * Clone the trajectory
     * @return a copy of this trajectory
     */
  //  virtual Trajectory<T, nx, nu, np, _ns>* Clone() const;
  
  /**
   * Updates the discrete trajectory using the parameter vector s, i.e. : s -> (ts,xs,us,p)
   * @return true if the resulting discrete trajectory is feasible
   */
  virtual bool Update(); 
  
  
  /**
   * Resize the discrete trajectory keeping any
   * previous states and if necessary adding empty states
   * at the end.
   * @param sn new number of discrete segments
   */ 
  //  void Resize(int sn);
  
  /**
   * Reverse the states+control alongs trajectory. keep same times
   */ 
  //    void Reverse();
  
  
  /**
   * Initialize the trajectory using sn segments
   * and optionally (if both si and sf are povided)
   * interpolate the states between
   * states si and sf. In fact, Init(sn) is equivalent to
   * Resize(sn).
   * @param sn number of segments
   * @param si start state (optional)
   * @param sf final state (optional)
   */
  // void Init(int sn,
  //    const State *si = 0,
  //    const State *sf = 0);
  
  
  /**
   * Attach trajectory traj to this trajectory. This trajectory
   * is modified and its size becomes as explained below. This operation
   * is fairly efficient since it is based on raw memory manipulation.
   * @param traj trajectory to be attached
   * @param back if true attach to back of this trajectory, else to the front
   * @param time if true then adjust time along the newly added states by
   *        incrementing/decrementing each state's time based on the 
   *        last/first state's time in the current trajectory depending
   *        on whether attaching to back/front, respectively.
   * @param js if true then treat the end point of this trajectory and 
   *           the first point of the given trajectory traj as the same 
   *           points (with same time)--in this case the total
   *           number of segments would be this->sn + traj.sn,
   *           otherwise it is this->sn + traj.sn + 1
  void Attach(const Trajectory<T, nx, nu, np, _ns> &traj,
              bool back = true,
              bool time = false,
              bool js = true);
   */
  
    /**
     * Append a state to the end of this trajectory
     * @param s state to be attached
     */
    //    void Append(const State &s);
  
  
  
    /**
     * Clear all states along this trajectory and free
     * associated memory. 
     */
  //  void Clear(); 
  
  /**
   * Set the path times starting from t0 with timestep h
   * 
   * @param t0 start time
   * @param h time step
   */
  //  void SetTime(double t0, double h);
  
  
  /**
   * Get the state at time t (t should already be set in the state s)
   * t should be in the valid range of trajectory times. This uses
   * the base interpolation define in System::Get. For now this method
   * assumes that all segments have equal time duration. 
   *
   * @param x state
   * @param u control
   * @param t time
   * @param p static parameter (optional)
   * @return true if arguments are feasible
   */
  //    void Get(State &s) const;
  virtual bool Get(T &x, Vectorcd &u, double t) const;
  
    
    /**
     * Get the index of trajectory segment in which time t falls. 
     * For now this method 
     * assumes that all segments have equal time duration.
     * @param t time
     * @return discrete trajectory segment where t falls (returns -1 if t
     * is out of bounds)
     */
    int Get(double t) const;


    /**
     * Get a subtrajectory given two starting times. Utilizes the System::Get
     * implementation.
     * @param traj to be filled in (the number of discrete segments traj.sn is used to determine the time step, so the discrete length of traj is not modified). Hence, one should pass traj with already existing states, and with traj->sn >0. If traj->sn=0 then just use GetState.
     * @param ti initial time (must be within this trajectory)
     * @param tf final time (must be within this trajectory)
     * @return true if arguments are feasible
     */
  // bool Get(Trajectory<T, nx, nu, np, _ns> &traj, double ti, double tf) const;
    
    
    /**
     * Add a point interpolated at s inside every trajectory segment
     * new trajectory has 2*sn+1 points
     * @param a number in [0,1] indicating where inside the segment to put the new point
     */
  // void Refine(double a = .5);
    
    /**
     * Resample the trajectory using sn new equaly spaced segments
     * @param sn new number of discrete segments
     */    
    // void Resample(int sn);


    /**
     * Add mn new states to the trajectory after points 
     * with indices in mi and interpolated between states at
     * indices mi[i] and mi[i+1] using the number mu[i] in [0,1].
     * @param mn number of new states to insert in the trajectory
     * @param mi (mn-array) indices of old states after which the new point will be inserted
     * @param mu (mn-array) numbers in the range [0,1] indicating where in the
     * mi[i]-th segment the new state will be interpolated
     */
  //    void Modify(int mn, const int *mi, const double *mu);


    /**
     * Checks whether this trajectory contains time t
     * @param t given time
     * @return true if t is between the start and end times of the trajectory
     */
  //    bool IsValidTime(double t) const;
    

    const System<T, nx, nu, np> &sys;    ///< mechanical system
  
  int ns; ///< trajectory parameter dimension

    Vectorsd s;           ///< continuous trajectory parametrization vector (this given vector should uniqely determine the whole continuous trajectory, typically together with a given initial state x0 and initial time t0). Examples of what s should be include: 1) simply a discrete sequence of constant controls applied during each discrete time-segment (this is the default -- in which case, s represents the sequence of controls us); 2) interpolating points for more general control curve parametrization using spline basis functions, 3) interpolating points for flat output parametrization using spline basis functions, 4) both control, state, and static params parametrizations often requires for direct collocation/multiple-shooting methods.

    vector<double> ts;    ///< sequence of times (length is N+1, and always at least 1 element, the first is the start time t0)
    vector<T> xs;        ///< sequence of states (length is N+1, and always at least 1, the first is the start state x0)
    vector<Vectorcd> us;  ///< sequence of controls (length is N), i.e. i-th control is regarded as the average control input over i-th segment
    Vectormd p;           ///< static parameters (could be none)
    
  //    bool ext;             ///< extended formulation (used for optimization purposes -- this flag is internally set when one requires variations in the exteneded stapce: space + time)

  //    double hmin;          ///< minimum timestep size (used internally for checking trajectory consistency) it is 1e-10 by default
    
    /**
     * Log the trajectory to a file.
     * @param logName log filename
     * @param di log every di-th state (optional, default is 1)
     */
  //    virtual void Log(const char *logName, int di = 1);

    /**
     * Log a state in simple human readable format
     * @param os output stream
     * @param di log every di-th state (optional, default is 1)
     */
  //    virtual void Log(std::ofstream& os, int di = 1);

    /**
     * Write a formatted number
     * @param os output stream
     * @param fw field width
     * @param a the number of write
     */
  //    static void Log(std::ofstream& os, int fw, double a);       


    
  private:
  //    friend std::ostream& operator<<(std::ostream &os, const Trajectory<T, nx, nu, np, _ns> &traj);
  //    friend std::istream& operator>>(std::istream &is, Trajectory &traj);
  };

  /**
   * Output the trajectory to a stream
   * @param os output stream
   * @param traj trajectory 
   * @return the output stream
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic,
    int _ns = Dynamic> 
  std::ostream& operator<<(std::ostream &os, const gcop::Trajectory<T, nx, nu, np, _ns> &traj);
   */


  /**
   * Input the trajectory from a stream
   * @param is input stream
   * @param traj trajectory
   * @return the input stream
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic,
    int _ns = Dynamic> 
  std::istream& operator>>(std::istream &is, gcop::Trajectory<T, nx, nu, np, _ns> &traj);
   */

  template <typename T, int nx, int nu, int np, int _ns> 
    Trajectory<T, nx, nu, np, _ntp>::Trajectory(System<T, nx, nu, np> &sys, int ns) : 
    sys(sys), ns(_ns != Dynamic ? _ns : ns), 
    ts(this->ns + 1), xs(this->ns + 1), us(this->ns), p(sys.P.n) {
    if (_ns == Dynamic)
      s.resize(this->ns);
  }
  

  template <typename T, int nx, int nu, int np, int _ns> 
    bool Trajectory<T, nx, nu, np, _ns>::Get(T &s, 
                                             const T &sa, 
                                             const T &sb, 
                                             double a) const
{
  Vectornd dx;
  sys.X.Lift(dx, sa, sb);
  sys.X.Retract(s, sa, a*dx);
}



  template <typename T, int nx, int nu, int np, int _ns> 
    void Trajectory<T, nx, nu, np, _ns>::Update() {
    assert(ns == us.size()*sys.U.n);
    for (int i = 0; i < us.size(); ++i) {
      memcpy(s.data() + i*sys.U.n, us[i].data(), sys.U.n*sizeof(double));
    }
    // assume parameter p does not play a role
  }
  
  template <typename T, int nx, int nu, int np, int _ns> 
    bool Trajectory<T, nx, nu, np, _ns>::Get(T &x, Vectorcd &u, double t) const {
    int N = us.size();
    if (!N)
      return;
    
    // find the time slot
    if ( (t < ts[0] - hmin) || (t > ts[N] + hmin)) {
      cerr << "[W] dcop::Trajectory::Get:\t time t= " << t << " out of bounds [" << ts[0] << "," << ts[N] << "]!" << endl;
      return;
    }
    
    
    // assume constant timestep  
    double h = ts[1] - ts[0];
    
    if (h < hmin) {
      cerr << "[W] dgc::Trajectory::GetState:\t time step is very small h=" << h << "!" << endl;    
    }
    
    
    if (fabs(t - ts[0]) < hmin) {
      x = xs[0];
      u = us[0];
      return;
    }
    
    if (fabs(t - ts[N]) < hmin) {
      x = xs[N];
      u = us[N-1];
      return;
    }


    double i = (t - ts[0])/h;
    int i0 = (int)i;
    
    //  cout << "t=" << s.t << " ti=" << states[0]->t << " tf=" << states[sn]->t << " h=" << h << endl;
    //  cout << " sn=" << sn << " i=" << i << " i0=" << i0 << endl;
    
    assert(i0 >= 0 && i0 < N);
    
    //  cout << "i0=" << i0 << " i=" << i <<  " sn=" << sn << endl;
    
    double a = i - i0;
    //    u = (1 - a)*us[i0] + a*us[i0 + 1]; 
    u = us[i0];
    Get(x, xs[i0], xs[i0+1], a);
  }
  
  
  template <typename T, int nx, int nu, int np, int _ns> 
    int Trajectory<T, nx, nu, np, _ns>::Get(double t) const {
    int N = us.size();
    // find the time slot
    if ( (t < ts[0] - hmin) || (t > ts[N] + hmin)) {
      cerr << "[W] gcop::Trajectory::Get:\t time t= " << t << " out of bounds [" << ts[0] << "," << ts[N] << "]!" << endl;
      return -1;
    }
    
    // assume constant timestep  
    double h = ts[1] - ts[0];
    
    if (h < hmin) {
      cerr << "[W] gcop::Trajectory::GetState:\t time step is very small h=" << h << "!" << endl;
    }
    
    if (fabs(t - ts[0]) < hmin)
      return 0;
    
    if (fabs(t - ts[N]) < hmin)
      return N-1;
    
    return (int)((t - ts[0])/h);    
  }
}

#endif
