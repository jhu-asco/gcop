#include "ros/ros.h"
#include <iomanip>
#include <iostream>
#include <dynamic_reconfigure/server.h>
#include "gcop_comm/CtrlTraj.h"//msg for publishing ctrl trajectory
#include <gcop/urdf_parser.h>
#include "tf/transform_datatypes.h"
#include <gcop/se3.h>
#include "gcop/dmoc.h" //gcop dmoc header
#include "gcop/lqcost.h" //gcop lqr header
#include "gcop/rn.h"
#include <tf/transform_listener.h>

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<MbsState> MbsDmoc;//defining chaindmoc

//ros messages
gcop_comm::CtrlTraj trajectory;

//Publisher
ros::Publisher trajpub;

//Timer
ros::Timer iteratetimer;

//Subscriber
//ros::Subscriber initialposn_sub;

//Pointer for mbs system
boost::shared_ptr<Mbs> mbsmodel;

//Pointer for Optimal Controller
boost::shared_ptr<MbsDmoc> mbsdmoc;

//MbsState final state
boost::shared_ptr<MbsState> xf;

int Nit = 1;//number of iterations for dmoc

void q2transform(geometry_msgs::Transform &transformmsg, Vector6d &bpose)
{
	tf::Quaternion q;
	q.setEulerZYX(bpose[2],bpose[1],bpose[0]);
	tf::Vector3 v(bpose[3],bpose[4],bpose[5]);
	tf::Transform tftransform(q,v);
	tf::transformTFToMsg(tftransform,transformmsg);
	//cout<<"---------"<<endl<<transformmsg.position.x<<endl<<transformmsg.position.y<<endl<<transformmsg.position.z<<endl<<endl<<"----------"<<endl;
}

void pubtraj() //N is the number of segments
{
	int N = mbsdmoc->us.size();
	cout<<"N: "<<N<<endl;
	int nb = mbsmodel->nb;
	cout<<"nb: "<<nb<<endl;
	Vector6d bpose;

	gcop::SE3::Instance().g2q(bpose, mbsdmoc->xs[0].gs[0]);
	q2transform(trajectory.statemsg[0].basepose,bpose);
	trajectory.statemsg[0].statevector.resize(nb-1);
	trajectory.statemsg[0].names.resize(nb-1);

	for(int count1 = 0;count1 < nb-1;count1++)
	{
		trajectory.statemsg[0].statevector[count1] = mbsdmoc->xs[0].r[count1];
		trajectory.statemsg[0].names[count1] = mbsmodel->joints[count1].name;
	}

	for (int i = 0; i < N; ++i) 
	{
		trajectory.statemsg[i+1].statevector.resize(nb-1);
		trajectory.statemsg[i+1].names.resize(nb-1);
		gcop::SE3::Instance().g2q(bpose, mbsdmoc->xs[i+1].gs[0]);
		q2transform(trajectory.statemsg[i+1].basepose,bpose);
		for(int count1 = 0;count1 < nb-1;count1++)
		{
			trajectory.statemsg[i+1].statevector[count1] = mbsdmoc->xs[i+1].r[count1];
			trajectory.statemsg[i+1].names[count1] = mbsmodel->joints[count1].name;
		}
		trajectory.ctrl[i].ctrlvec.resize(6+nb-1);
		for(int count1 = 0;count1 < 6+nb-1;count1++)
		{
			trajectory.ctrl[i].ctrlvec[count1] = mbsdmoc->us[i](count1);
		}
	}
	//final goal:
	gcop::SE3::Instance().g2q(bpose, xf->gs[0]);
	q2transform(trajectory.finalgoal.basepose,bpose);
	trajectory.finalgoal.statevector.resize(nb-1);
	trajectory.finalgoal.names.resize(nb-1);

	for(int count1 = 0;count1 < nb-1;count1++)
	{
		trajectory.finalgoal.statevector[count1] = xf->r[count1];
		trajectory.finalgoal.names[count1] = mbsmodel->joints[count1].name;
	}

	trajpub.publish(trajectory);

}
void iterateCallback(const ros::TimerEvent & event)
{
	//ros::Time startime = ros::Time::now(); 
	for (int count = 1;count <= Nit;count++)
		mbsdmoc->Iterate();//Updates us and xs after one iteration
	//double te = 1e6*(ros::Time::now() - startime).toSec();
	//cout << "Time taken " << te << " us." << endl;
	//publish the message
	pubtraj();
}
/*
void initialposnCallback(const geometry_msgs::TransformStamped::ConstPtr &currframe)
{
	tf::StampedTransform UV_O;
	transformStampedMsgToTF(*currframe,UV_O);//converts to the right format 
	//getrpy:

	double roll,pitch,yaw;
	UV_O.getBasis().getRPY(roll,pitch,yaw);
	double tcurr = currframe->header.stamp.toSec();
	tf::Vector3 y = UV_O.getOrigin();
	Vector4d x0 = Vector4d::Zero();// initial state
	x0[0] = y[0];
	x0[1] = y[1];
	x0[2] = yaw;
	x0[3] =0;/// need to calculate velocity
	xs[0] = x0;
	//ros::TimerEvent e1;
	//iterateCallback(e1);
	return;
}
*/

int main(int argc, char** argv)
{
	ros::init(argc, argv, "chainload");
	ros::NodeHandle n;
	//Initialize publisher
	trajpub = n.advertise<gcop_comm::CtrlTraj>("ctrltraj",2);
	//Subscribe to initial posn from tf
	//initialposn_sub = rosdmoc.subscribe("mocap",1,initialposnCallback);
	//get parameter for xml_string:
	string xml_string, xml_filename;
	if(!ros::param::get("/robot_description", xml_string))
	{
		ROS_ERROR("Could not fetch xml file name");
		return 0;
	}
	//Create Mbs system
	mbsmodel = gcop_urdf::mbsgenerator(xml_string);
	//set no gravity:
	mbsmodel->ag << 0, 0, -9.81;

	//Printing the mbsmodel params:
	for(int count = 0;count<(mbsmodel->nb);count++)
	{
		cout<<"Ds["<<mbsmodel->links[count].name<<"]"<<endl<<mbsmodel->links[count].ds<<endl;
		cout<<"I["<<mbsmodel->links[count].name<<"]"<<endl<<mbsmodel->links[count].I<<endl;
	}
	for(int count = 0;count<(mbsmodel->nb)-1;count++)
	{
		cout<<"Joint["<<mbsmodel->joints[count].name<<"].gc"<<endl<<mbsmodel->joints[count].gc<<endl;
		cout<<"Joint["<<mbsmodel->joints[count].name<<"].gp"<<endl<<mbsmodel->joints[count].gp<<endl;
		cout<<"Joint["<<mbsmodel->joints[count].name<<"].a"<<endl<<mbsmodel->joints[count].a<<endl;
	}

	//Using it:
	//define parameters for the system
	int nb = mbsmodel->nb;
	int N = 100;      // discrete trajectory segments
	double tf = 20;   // time-horizon
	double h = tf/N; // time-step

	//times
	vector<double> ts(N+1);
	for (int k = 0; k <=N; ++k)
		ts[k] = k*h;


	//Define Final State
	xf.reset( new MbsState(nb));
	xf->gs[0].setIdentity();
	//xf->gs[0](0,3) = 0.0;
	xf->vs[0].setZero();
	xf->dr.setZero();
	xf->r[0] =0;
	xf->r[1] =0;
	//xf->r.fill(M_PI/4);   

	//Define Lqr Cost
	LqCost<MbsState> cost(mbsmodel->X, (Rn<>&)mbsmodel->U, tf, *xf);
	cost.Qf(0,0) = 2; cost.Qf(1,1) = 2; cost.Qf(2,2) = 2;
	cost.Qf(3,3) = 50; cost.Qf(4,4) = 50; cost.Qf(5,5) = 50;
	//cost.Qf(9,9) = 20; cost.Qf(10,10) = 20; cost.Qf(11,11) = 20;
	
	
//Initial State
	MbsState x(nb);
  x.gs[0].setIdentity();
	//x.gs[0](0,3) = 1.0;
	x.r[0] = -M_PI/2;
	x.r[1] = M_PI/2;
  //x.r.setZero();   
  mbsmodel->FK(x);

  // states
  vector<MbsState> xs(N+1, x);
  xs[0].vs[0].setZero();
//	xs[0].r[0] = 1.57;
	//xs[0].r[1] = -1.57;
  xs[0].dr[0] = 0;
  xs[0].dr[1] = 0;
	//cout<<"xs0"<<xs[0].r[0]<<"\t"<<xs[0].r[1]<<"\t"<<endl;
	mbsmodel->KStep(xs[0], x, h);




	//Define the initial state mbs
	/*
	MbsState x(nb);
	x.gs[0].setIdentity();
	x.gs[0](0,3) = -2;
	//gcop::SE3::Instance().rpyxyz2g(x.gs[0], Vector3d(0,0,0), Vector3d(-2,0,0));
	x.r[0] = -M_PI/2;
	x.r[1] = M_PI/2;
	x.vs[0].setZero();
	x.dr.fill(0);
	mbsmodel->Rec(x, h);
	cout<<"x.r0: "<<x.r[0]<<"\t"<<x.r[1]<<endl;
	*/


	// initial controls (e.g. hover at one place)
	VectorXd u(6+ (nb-1));
	u.setZero();
	u[5] = 30*9.81;

	//States and controls for system
	vector<VectorXd> us(N,u);
	//vector<MbsState> xs(N+1,x);

	mbsdmoc.reset(new MbsDmoc(*mbsmodel, cost, ts, xs, us));
	mbsdmoc->mu = 10;

	//Trajectory message initialization
	trajectory.N = N;
	trajectory.statemsg.resize(N+1);
	trajectory.ctrl.resize(N);
	trajectory.time = ts;
	trajectory.finalgoal.statevector.resize(nb-1);
	
	//Debug true for mbs

  mbsmodel->debug = true;

	// Create timer for iterating	
	iteratetimer = n.createTimer(ros::Duration(0.02), iterateCallback);
	iteratetimer.start();
	ros::spin();
	return 0;
}


	




