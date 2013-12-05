#include "ros/ros.h"
#include <iomanip>
#include <iostream>
#include "gcop/dmoc.h" //gcop dmoc header
#include "gcop/rnlqcost.h" //gcop lqr header
#include "gcop/rccar.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<Vector4d, 4, 2> RccarDmoc;


int main(int argc, char** argv)
{
	ros::init(argc, argv, "rccarctrl");
	ros::NodeHandle mavlink_nh("/rcctrl");

	//define parameters
  int N = 64;        // number of segments
  double tf = 20,mu = 0.01;    // time horizon
  Vector4d x0 = Vector4d::Zero();// initial state
  Vector4d xf = Vector4d::Zero();// final state
  VectorXd Q(4);//Costs
  VectorXd R(2);  
  VectorXd Qf(4);

	//get parameters from ros:
	ros::param::get("/dmoc/tf", tf);
	ros::param::get("/dmoc/N", N);

	ros::param::get("/dmoc/x0", x0(0));
	ros::param::get("/dmoc/y0", x0(1));
	ros::param::get("/dmoc/vx0", x0(2));
	ros::param::get("/dmoc/vy0", x0(3));

	ros::param::get("/dmoc/xN", xf(0));
	ros::param::get("/dmoc/xN", xf(1));
	ros::param::get("/dmoc/vy0", xf(2));
	ros::param::get("/dmoc/vyN", xf(3));

	ros::param::get("/dmoc/Qf1", Qf(0));
	ros::param::get("/dmoc/Qf2", Qf(1));
	ros::param::get("/dmoc/Qf3", Qf(2));
	ros::param::get("/dmoc/Qf4", Qf(3));

	ros::param::get("/dmoc/Q1", Q(0));
	ros::param::get("/dmoc/Q2", Q(1));
	ros::param::get("/dmoc/Q3", Q(2));
	ros::param::get("/dmoc/Q4", Q(3));

	ros::param::get("/dmoc/R1", R(0));
	ros::param::get("/dmoc/R2", R(1));

	ros::param::get("/dmoc/mu", mu);

  //conversions:
  double h = tf/N;   // time step
  // cost
  RnLqCost<4, 2> cost(tf, xf);

	cost.Q = Q.asDiagonal();
  cost.R = R.asDiagonal();
  cost.Qf = Qf.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Vector4d> xs(N+1);
  // initial state
  xs[0] = x0;

  // initial controls
  vector<Vector2d> us(N);
  for (int i = 0; i < N/2; ++i) {
    us[i] = Vector2d(.01, .0);
    us[N/2+i] = Vector2d(-.01, .0);
  }
 
 
  //Initializing systems  
  Rccar sys;
  RccarDmoc dmoc(sys, cost, ts, xs, us);  
  dmoc.mu = mu;

  //struct timeval timer;
  // dmoc.debug = false; // turn off debug for speed
  getchar();

  for (int i = 0; i < 30; ++i) {
    ros::Time startime = ros::Time::now(); 
    dmoc.Iterate();
		double te = 1e6*(ros::Time::now() - startime).toSec();
    cout << "Iteration #" << i << " took: " << te << " us." << endl;
    getchar();
  }

  cout << "done!" << endl;
  while(1)


  return 0;
}
