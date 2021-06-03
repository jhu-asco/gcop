
#ifndef GCOP_BODY3DWAYPOINTCOST_H
#define GCOP_BODY3DWAYPOINTCOST_H

#include "body3d.h"
#include "body3dcost.h"
#include "lqcost.h"
#include "multicost.h"
#include <iostream>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  /**
   * Linear-quadratic cost for rigid body systems
   *
   * Author: Marin Kobilarov
   */
   class Body3dWaypointCost : public Cost<Body3dState,12,4> {

    typedef Matrix<double, 12, 12> Matrixnd;
    typedef Matrix<double, 4, 4> Matrixcd;

  public:
    
    Body3dWaypointCost(Body3d<4> &sys, vector<double> &tfs, vector<Body3dState> &goals);
    virtual double L(double t, const Body3dState& x, const Vector4d& u, double h,
                     const Matrix<double, Dynamic, 1> *p = 0,
                     Matrix<double, 12, 1> *Lx = 0, Matrix<double, 12, 12>* Lxx = 0,
                     Matrix<double, 4, 1> *Lu = 0, Matrix<double, 4, 4>* Luu = 0,
                     Matrix<double, 12, 4> *Lxu = 0,
                     Matrix<double,Dynamic, 1> *Lp = 0, Matrix<double,Dynamic,Dynamic> *Lpp = 0,
                     Matrix<double, Dynamic, 12> *Lpx = 0);

    //vector<Body3dCost<4>> cost_list;
    vector<Body3dState> goal_list;
    vector<double> time_list;
    Matrixnd Q;
    Body3d<4> &m_system;
    //Body3dCost<4> m_cost;
  };  
  
  Body3dWaypointCost::Body3dWaypointCost(Body3d<4> &sys, vector<double> &tfs, vector<Body3dState> &goals): 
    Cost<Body3dState,12,4>(sys,tfs.back()),m_system(sys) {

    //cout << "Initialization Entered" << endl;
    //vector<double> time_list;
    int num_waypoints = goals.size();
    //cout << "Num _waypoints: " << num_waypoints << endl;
    if (tfs.size() != num_waypoints) {
        double tf = tfs.back();
	//cout << "Tf: " << tf << endl;
        double dt = tf/num_waypoints;
	//cout << "dt: " << dt << endl;
        time_list = vector<double>(num_waypoints);
        for (int ii = 0; ii < num_waypoints; ++ii) {
            time_list[ii] = dt*(ii+1);
	    //cout << "time list [" << ii << "] set to " << time_list[ii] << endl;
        }
    } else {
        time_list = tfs;
    }
    //cout << "time_list length now is " << time_list.size() << endl;
    //double tf = time_list.back();
    //Body3dCost<4> temp(sys,tf,goals.back()
    /*Body3dState last_goal = goals.back();
    m_cost(sys,tf,last_goal);
    vector<Body3dState> xds(N+1);
    vector<Vector4d> uds(N);
    double dt = tf/N;
    xds[0] = goals[0];
    for(int ii = 0; ii < N; ++ii){
	    uds[ii].head(3).setZero();
	    uds[ii][3] = 9.81*sys.m;
	    double time = dt*(ii+1);
	    for(int jj = 0; jj < num_waypoints; ++jj){
		    if (time < (time_list[jj]-1e-10)) {
			    xds[ii+1] = goals[jj];
			    break;
		    }
		    xds[ii+1] = goals.back();
	    }
    }
    m_cost.SetReference(&xds,&uds);*/
    goal_list = goals;
    Q.setIdentity();
    /*double previous_time = 0;
    for(int ii = 0; ii < num_waypoints; ++ii){
        double dt = time_list[ii] - previous_time;
        Body3dState goal = goals[ii];
        Body3dCost<4> temp_cost(sys,dt,goal);
        cost_list.push_back(temp_cost);
	previous_time = time_list[ii];
    }*/
    //cout << "Initialization Finished" << endl;
  }
  double Body3dWaypointCost::L(double t, const Body3dState& x, const Vector4d& u, double h,
                                          const Matrix<double, Dynamic, 1> *p,
                                          Matrix<double, 12, 1> *Lx, Matrix<double, 12, 12>* Lxx,
                                          Matrix<double, 4, 1> *Lu, Matrix<double, 4, 4>* Luu,
                                          Matrix<double, 12, 4> *Lxu,
                                          Matrix<double, Dynamic, 1> *Lp, Matrix<double,Dynamic,Dynamic> *Lpp,
                                          Matrix<double, Dynamic, 12> *Lpx) {
    /*double previous_time = 0;
    for(int ii = 0; ii < cost_list.size(); ++ii) {
        if (t < time_list[ii]) {
            return cost_list[ii].L(t-previous_time, x, u, h, p, Lx, Lxx, Lu, Luu, Lxu, Lp, Lpp, Lpx);
        }
        previous_time = time_list[ii];
    }
    //If we reach here, we've reached the end of the time_list, so use the last one
    int num_waypoints = time_list.size();
    Body3dCost<4> last_cost = cost_list.back();
    double next_to_last_time = time_list[num_waypoints-2];
    return last_cost.L(t-next_to_last_time, x, u, h, p, Lx, Lxx, Lu, Luu, Lxu, Lp, Lpp, Lpx);
    */
    for(int ii = 0; ii < goal_list.size(); ++ii) {
	    double wp_time = time_list[ii];
	    //cout << "Iterating: " << ii << " wp time: " << wp_time << " t: " << t<< endl;
	    if ((t > wp_time - 1e-5) && (t < wp_time + 1e-5)) {
		    //cout << "Time " << t << "Triggered wp time " << wp_time << endl;
		    Body3dState goal = goal_list[ii];
		    Body3dCost<4> temp(m_system,t+1e-11,goal);
		    temp.Qf = Q;
                    double temp_L = temp.L(t,x,u,h,p,Lx,Lxx,Lu,Luu,Lxu,Lp,Lpp,Lpx);
                    if (std::isnan(temp_L)) {
                      cout << "NaN t: " << t << endl;
                      cout << "WaypointCost threw NaN" << endl;
                      cout << "x: " << endl << "R: " << endl << x.R << endl << "p: " << endl << x.p << endl << "w: " << endl << x.w << endl << "v: " << endl << x.v << endl;
                      cout << "xf: " << endl << "R: " << endl << goal.R << endl << "p: " << endl << goal.p << endl << "w: " << endl << goal.w << endl << "v: " << endl << goal.v << endl;
                    }
		    return temp_L;
	    }
    }
    //If here, we reached the end of time_list
    if (Lx)
      Lx->setZero();
    if (Lxx)
      Lxx->setZero();
    if (Lu)
      Lu->setZero();
    if (Luu)
      Luu->setZero();
    if (Lxu)
      Lxu->setZero();
    if (Lp)
      Lp->setZero();
    if (Lpp)
      Lpp->setZero();
    if (Lpx)
      Lpx->setZero();
    return 0;

    //return m_cost.L(t,x,u,h,p,Lx,Lxx,Lu,Luu,Lxu,Lp,Lpp,Lpx);
  }
}


#endif
