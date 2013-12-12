#include "ros/ros.h"
#include <iomanip>
#include <iostream>
#include <dynamic_reconfigure/server.h>
#include "gcop/dmoc.h" //gcop dmoc header
#include "gcop/rnlqcost.h" //gcop lqr header
#include "gcop/rccar.h"
#include "gcop_comm/CtrlTraj.h"//msg for publishing ctrl trajectory
#include "gcop_ctrl/DMocInterfaceConfig.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<Vector4d, 4, 2> RccarDmoc;

//ros messages
gcop_comm::CtrlTraj trajectory;

//Publisher
ros::Publisher trajpub;

//Timer
ros::Timer iteratetimer;


//Initialize the rccar system
Rccar sys;

//Optimal Controller
	RccarDmoc *dmoc;

//Cost class
RnLqCost<4, 2>*cost;

//Define states and controls for system
	vector<double> ts;
  vector<Vector4d> xs;
  vector<Vector2d> us;
  Vector4d xf = Vector4d::Zero();// final state initialization passed by reference so needs to be global to be changed
	int Nit = 10;//number of iterations for dmoc


void pubtraj() //N is the number of segments
{
	int N = us.size();
	
	for (int count = 0;count<N+1;count++)
	{
		for(int count1 = 0;count1 < 4;count1++)
		{
			trajectory.statemsg[count].statevector.resize(4);

			trajectory.statemsg[count].statevector[count1] = xs[count](count1);
		}
	}
	for (int count = 0;count<N;count++)
	{
		for(int count1 = 0;count1 < 2;count1++)
		{
			trajectory.ctrl[count].ctrlvec.resize(2);

			trajectory.ctrl[count].ctrlvec[count1] = us[count](count1);
		}
	}
	trajectory.time = ts;
	//final goal:
	for(int count = 0;count<4;count++)
	{
		trajectory.finalgoal.statevector[count] = xf(count);
	}
	trajpub.publish(trajectory);
}

void iterateCallback(const ros::TimerEvent & event)
{
	//ros::Time startime = ros::Time::now(); 
	for (int count = 1;count <= Nit;count++)
		dmoc->Iterate();//Updates us and xs after one iteration
	//double te = 1e6*(ros::Time::now() - startime).toSec();
	//cout << "Time taken " << te << " us." << endl;
	//publish the message
	pubtraj();
}


void paramreqcallback(gcop_ctrl::DMocInterfaceConfig &config, uint32_t level) 
{
	Nit = config.Nit; 
	int N = config.N;
	double h = config.tf/N;   // time step
	//Vector4d xf = Vector4d::Zero();// final state
	Vector4d x0 = Vector4d::Zero();// initial state
	VectorXd Q(4);//Costs
	VectorXd R(2);  
	VectorXd Qf(4);

	xf(0) = config.xN;
	xf(1) = config.yN;
	xf(2) = config.thetaN;
	xf(3) = config.vN;

	x0(0) = config.x0;
	x0(1) = config.y0;
	x0(2) = config.theta0;
	x0(3) = config.v0;

	Q(0) = config.Q1;
	Q(1) = config.Q2;
	Q(2) = config.Q3;
	Q(3) = config.Q4;

	Qf(0) = config.Qf1;
	Qf(1) = config.Qf2;
	Qf(2) = config.Qf3;
	Qf(3) = config.Qf4;

	R(0) = config.R1;
	R(1) = config.R2;

	//resize
	ts.resize(N+1);
	xs.resize(N+1);
	us.resize(N);
	trajectory.N = N;
	trajectory.statemsg.resize(N+1);
	trajectory.ctrl.resize(N);
	
	// cost

	//RnLqCost<4, 2> cost(config.tf, xf);
	cost->tf = config.tf;
	//cost->xf = xf; //xf is already changed when we assign xf its new values

	cost->Q = Q.asDiagonal();
	cost->R = R.asDiagonal();
	cost->Qf = Qf.asDiagonal();


	for (int k = 0; k <=N; ++k)
		ts[k] = k*h;

	// initial state
	xs[0] = x0;

	//initial controls
	/*for (int i = 0; i < N/2; ++i) {
		us[i] = Vector2d(.01, .0);
		us[N/2+i] = Vector2d(-.01, .0);
	}
	*/
	//change parameters in dmoc:
	dmoc->ts = ts;
//	dmoc->xs = xs;
	//dmoc->us = us;
	dmoc->mu = config.mu;



	//dont know what to do with the cost


	//destroy previous dmoc:
	//delete(dmoc);
	//dmoc = new RccarDmoc(sys, cost, ts, xs, us);  
}

int main(int argc, char** argv)
{
	ros::init(argc, argv, "rccarctrl");
	ros::NodeHandle rosdmoc("/dmoc");
	//Initialize publisher
	trajpub = rosdmoc.advertise<gcop_comm::CtrlTraj>("ctrltraj",1);

	//define parameters for the system
  int N = 64;        // number of segments
  double tf = 20,mu = 0.01;    // time horizon
  Vector4d x0 = Vector4d::Zero();// initial state
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

	//resize the states and controls
	ts.resize(N+1);
	xs.resize(N+1);
	us.resize(N);

	  //conversions:
  double h = tf/N;   // time step

	cost = new RnLqCost<4, 2>(tf,xf);

	cost->Q = Q.asDiagonal();
  cost->R = R.asDiagonal();
  cost->Qf = Qf.asDiagonal();

  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // initial state
  xs[0] = x0;
  //initial controls
  for (int i = 0; i < N/2; ++i) {
    us[i] = Vector2d(.01, .0);
    us[N/2+i] = Vector2d(-.01, .0);
  }
 
 
  dmoc = new RccarDmoc(sys, *cost, ts, xs, us);  
  dmoc->mu = mu;


  //Trajectory message initialization
	trajectory.N = N;
	trajectory.statemsg.resize(N+1);
	trajectory.ctrl.resize(N);
	trajectory.time = ts;
	trajectory.finalgoal.statevector.resize(4);
	//trajectory.time.resize(N);
//Dynamic Reconfigure setup Callback ! immediately gets called with default values	
	dynamic_reconfigure::Server<gcop_ctrl::DMocInterfaceConfig> server;
	dynamic_reconfigure::Server<gcop_ctrl::DMocInterfaceConfig>::CallbackType f;
	f = boost::bind(&paramreqcallback, _1, _2);
	server.setCallback(f);
	//create timer for iteration
  iteratetimer = rosdmoc.createTimer(ros::Duration(0.02), iterateCallback);
	iteratetimer.start();
	ros::spin();
  return 0;
}
