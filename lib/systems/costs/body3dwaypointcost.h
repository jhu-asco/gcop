
#ifndef GCOP_BODY3DWAYPOINTCOST_H
#define GCOP_BODY3DWAYPOINTCOST_H

#include "body3d.h"
#include "body3dcost.h"
#include "lqcost.h"
#include "multicost.h"

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

    vector<Body3dCost<4>*> cost_list;
    vector<double> time_list;
  };  
  
  Body3dWaypointCost::Body3dWaypointCost(Body3d<4> &sys, vector<double> &tfs, vector<Body3dState> &goals): 
    Cost<Body3dState,12,4>(sys,tfs.back()) {

    int num_waypoints = goals.size();
    if (tfs.size() != num_waypoints) {
        double tf = tfs.back();
        double dt = tf/num_waypoints;
        time_list = vector<double>(num_waypoints);
        for (int ii = 0; ii < num_waypoints; ++ii) {
            time_list[ii] = dt*(ii+1);
        }
    } else {
        time_list = tfs;
    }

    double previous_time = 0;
    for(int ii = 0; ii < num_waypoints; ++ii){
        double dt = time_list[ii] - previous_time;
        Body3dState goal = goals[ii];
        Body3dCost<4> temp_cost(sys,dt,goal);
        cost_list.push_back(&temp_cost);
    }
  }
  double Body3dWaypointCost::L(double t, const Body3dState& x, const Vector4d& u, double h,
                                          const Matrix<double, Dynamic, 1> *p,
                                          Matrix<double, 12, 1> *Lx, Matrix<double, 12, 12>* Lxx,
                                          Matrix<double, 4, 1> *Lu, Matrix<double, 4, 4>* Luu,
                                          Matrix<double, 12, 4> *Lxu,
                                          Matrix<double, Dynamic, 1> *Lp, Matrix<double,Dynamic,Dynamic> *Lpp,
                                          Matrix<double, Dynamic, 12> *Lpx) {
    double previous_time = 0;
    for(int ii = 0; ii < cost_list.size(); ++ii) {
        if (t < time_list[ii]) {
            return cost_list[ii]->L(t-previous_time, x, u, h, p, Lx, Lxx, Lu, Luu, Lxu, Lp, Lpp, Lpx);
        }
        previous_time = time_list[ii];
    }
    //If we reach here, we've reached the end of the time_list, so use the last one
    int num_waypoints = time_list.size();
    return cost_list.back()->L(t-time_list[num_waypoints-2], x, u, h, p, Lx, Lxx, Lu, Luu, Lxu, Lp, Lpp, Lpx);
  }
}


#endif
