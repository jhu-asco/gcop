#include "ros/ros.h"
#include <iomanip>
#include <iostream>
#include <dynamic_reconfigure/server.h>
#include "gcop_comm/CtrlTraj.h"//msg for publishing ctrl trajectory
#include <gcop/urdf_parser.h>
#include "tf/transform_datatypes.h"
#include <gcop/se3.h>
#include "gcop/mbscontroller.h"
#include "gcop_ctrl/MbsSimInterfaceConfig.h"
#include <tf/transform_listener.h>
#include <XmlRpcValue.h>

using namespace std;
using namespace Eigen;
using namespace gcop;



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


//States and controls for system
vector<VectorXd> us;
vector<MbsState> xs;
vector<double> ts;

string mbstype; // Type of system
double tfinal = 20;   // time-horizon
Matrix4d gposeroot_i; //inital inertial frame wrt the joint frame

void rpy2transform(geometry_msgs::Transform &transformmsg, Vector6d &bpose)
{
	tf::Quaternion q;
	q.setEulerZYX(bpose[2],bpose[1],bpose[0]);
	tf::Vector3 v(bpose[3],bpose[4],bpose[5]);
	tf::Transform tftransform(q,v);
	tf::transformTFToMsg(tftransform,transformmsg);
	//cout<<"---------"<<endl<<transformmsg.position.x<<endl<<transformmsg.position.y<<endl<<transformmsg.position.z<<endl<<endl<<"----------"<<endl;
}

void xml2vec(VectorXd &vec, XmlRpc::XmlRpcValue &my_list)
{
	ROS_ASSERT(my_list.getType() == XmlRpc::XmlRpcValue::TypeArray);
	ROS_ASSERT(my_list.size() > 0);
	vec.resize(my_list.size());
	//cout<<my_list.size()<<endl;
	//ROS_ASSERT(vec.size() <= my_list.size()); //Desired size

	for (int32_t i = 0; i < my_list.size(); i++) 
	{
				ROS_ASSERT(my_list[i].getType() == XmlRpc::XmlRpcValue::TypeDouble);
				cout<<"my_list["<<i<<"]\t"<<my_list[i]<<endl;
			  vec[i] =  (double)(my_list[i]);
	}
}
void simtraj(const ros::TimerEvent &event) //N is the number of segments
{
	int N = us.size();
	double h = tfinal/N; // time-step
	//cout<<"N: "<<N<<endl;
	int csize = mbsmodel->U.n;
	//cout<<"csize: "<<csize<<endl;
	int nb = mbsmodel->nb;
	//cout<<"nb: "<<nb<<endl;
	Vector6d bpose;

	gcop::SE3::Instance().g2q(bpose,gposeroot_i*xs[0].gs[0]);
	rpy2transform(trajectory.statemsg[0].basepose,bpose);

	for(int count1 = 0;count1 < nb-1;count1++)
	{
		trajectory.statemsg[0].statevector[count1] = xs[0].r[count1];
		trajectory.statemsg[0].names[count1] = mbsmodel->joints[count1].name;
	}

	for (int i = 0; i < N; ++i) 
	{
//		cout<<"i "<<i<<endl;
		mbsmodel->Step(xs[i+1], i*h, xs[i], us[i], h);
		//getchar();
		gcop::SE3::Instance().g2q(bpose, gposeroot_i*xs[i+1].gs[0]);
		rpy2transform(trajectory.statemsg[i+1].basepose,bpose);
		for(int count1 = 0;count1 < nb-1;count1++)
		{
			trajectory.statemsg[i+1].statevector[count1] = xs[i+1].r[count1];
			trajectory.statemsg[i+1].names[count1] = mbsmodel->joints[count1].name;
		}
		for(int count1 = 0;count1 < csize;count1++)
		{
			trajectory.ctrl[i].ctrlvec[count1] = us[i](count1);
		}
	}


	trajectory.time = ts;

	trajpub.publish(trajectory);

}


void paramreqcallback(gcop_ctrl::MbsSimInterfaceConfig &config, uint32_t level) 
{
	int nb = mbsmodel->nb;

	if(level == 0xffffffff)
		config.tf = tfinal;

	tfinal = config.tf;

	if(level & 0x00000001)
	{

		if(config.i_J > nb-1)
			config.i_J = nb-1; 
		else if(config.i_J < 1)
			config.i_J = 1;

		config.Ji = xs[0].r[config.i_J-1];     
		config.Jvi = xs[0].dr[config.i_J-1];     

		if(config.i_u > mbsmodel->U.n)
			config.i_u = mbsmodel->U.n;
		else if(config.i_u < 1)
			config.i_u = 1;
		config.ui = us[0][config.i_u-1];
	}

/*	if(level & 0xffffffff)
		config.N = N;

	if(level & 0x00000004)
	{
		cout<<"I am called"<<endl;
		//resize
		N = config.N;
		ts.resize(N+1);
		xs.resize(N+1,xs[0]);
		us.resize(N,us[0]);
		trajectory.N = N;
		trajectory.statemsg.resize(N+1);
		trajectory.ctrl.resize(N);

		for (int i = 0; i < N; ++i) 
		{
			trajectory.statemsg[i+1].statevector.resize(nb-1);
			trajectory.statemsg[i+1].names.resize(nb-1);
			trajectory.ctrl[i].ctrlvec.resize(mbsmodel->U.n);
		}
	}
	*/

	int N = us.size();
	double h = tfinal/N;
	//Setting Values
	for (int k = 0; k <=N; ++k)
		ts[k] = k*h;

  if((mbstype == "FLOATBASE")||(mbstype == "AIRBASE"))
	{
		gcop::SE3::Instance().rpyxyz2g(xs[0].gs[0], Vector3d(config.roll,config.pitch,config.yaw), Vector3d(config.x,config.y,config.z));
		xs[0].vs[0]<<config.vroll, config.vpitch, config.vyaw, config.vx, config.vy, config.vz;
	}

	xs[0].r[config.i_J -1] = config.Ji;
	xs[0].dr[config.i_J -1] = config.Jvi;
	mbsmodel->Rec(xs[0], h);

	for(int count = 0;count <us.size();count++)
		us[count][config.i_u-1] = config.ui; 

	return;
}


int main(int argc, char** argv)
{
	ros::init(argc, argv, "chainload");
	ros::NodeHandle n("mbsdmoc");
	//Initialize publisher
	trajpub = n.advertise<gcop_comm::CtrlTraj>("ctrltraj",1);
	//get parameter for xml_string:
	string xml_string, xml_filename;
	if(!ros::param::get("/robot_description", xml_string))
	{
		ROS_ERROR("Could not fetch xml file name");
		return 0;
	}
	n.getParam("basetype",mbstype);
	assert((mbstype == "FIXEDBASE")||(mbstype == "AIRBASE")||(mbstype == "FLOATBASE"));
	VectorXd xmlconversion;
	Matrix4d gposei_root; //inital inertial frame wrt the joint frame
	//Create Mbs system
	mbsmodel = gcop_urdf::mbsgenerator(xml_string,gposei_root, mbstype);
	cout<<"gposei_root: "<<endl<<gposei_root<<endl;
	gcop::SE3::Instance().inv(gposeroot_i,gposei_root);
	//getchar();
	//mbsmodel->debug = true;
	//mbsmodel->method = 1;
	//mbsmodel->iters = 4;
	cout<<"Mbstype: "<<mbstype<<endl;
	mbsmodel->ag << 0, 0, -0.05;//default
	//get ag from parameters
	XmlRpc::XmlRpcValue ag_list;
	if(n.getParam("ag", ag_list))
		xml2vec(xmlconversion,ag_list);
	ROS_ASSERT(xmlconversion.size() == 3);
	//mbsmodel->ag = gposei_root.topLeftCorner(3,3).transpose()*xmlconversion.head(3);
	mbsmodel->ag = xmlconversion.head(3);

	cout<<"mbsmodel->ag: "<<endl<<mbsmodel->ag<<endl;
	//mbsmodel->ag = xmlconversion.head(3);

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

	cout<<"mbsmodel damping: "<<mbsmodel->damping.transpose()<<endl;
	int ctrlveclength = (mbstype == "AIRBASE")? 4:
		(mbstype == "FLOATBASE")? 6:
		(mbstype == "FIXEDBASE")? 0:0;
	cout<<"ctrlveclength "<<ctrlveclength<<endl;
	//set bounds on basebody controls
	XmlRpc::XmlRpcValue ub_list;
	if(n.getParam("ulb", ub_list))
	{
		xml2vec(xmlconversion,ub_list);
		ROS_ASSERT(xmlconversion.size() == ctrlveclength);
		mbsmodel->U.lb.head(ctrlveclength) = xmlconversion;
	}
	//upper bound
	if(n.getParam("uub", ub_list))
	{
		xml2vec(xmlconversion,ub_list);
		ROS_ASSERT(xmlconversion.size() == ctrlveclength);
		mbsmodel->U.ub.head(ctrlveclength) = xmlconversion;
	}
	cout<<"mbsmodel U.lb: "<<mbsmodel->U.lb.transpose()<<endl;
	cout<<"mbsmodel U.ub: "<<mbsmodel->U.ub.transpose()<<endl;
	cout<<"mbsmodel X.lb.gs[0]"<<endl<<mbsmodel->X.lb.gs[0]<<endl;
	cout<<"mbsmodel X.ub.gs[0]"<<endl<<mbsmodel->X.ub.gs[0]<<endl;
	cout<<"mbsmodel X.lb.vs[0]: "<<mbsmodel->X.lb.vs[0].transpose()<<endl;
	cout<<"mbsmodel X.ub.vs[0]: "<<mbsmodel->X.ub.vs[0].transpose()<<endl;
	//Initialize after setting everything up
	mbsmodel->Init();

	//Using it:
	//define parameters for the system
	int nb = mbsmodel->nb;

	//times
	int N = 100;      // discrete trajectory segments
	n.getParam("tf",tfinal);
	n.getParam("N",N);
	double h = tfinal/N; // time-step



	//Define the initial state mbs
	MbsState x(nb,(mbstype == "FIXEDBASE"));
	x.gs[0].setIdentity();
	x.vs[0].setZero();
	x.dr.setZero();
	x.r.setZero();

	// Get X0	 from params
	XmlRpc::XmlRpcValue x0_list;
	if(n.getParam("X0", x0_list))
	{
		xml2vec(xmlconversion,x0_list);
		ROS_ASSERT(xmlconversion.size() == 12);
		x.vs[0] = xmlconversion.tail<6>();
		gcop::SE3::Instance().rpyxyz2g(x.gs[0],xmlconversion.head<3>(),xmlconversion.segment<3>(3)); 
	}
	x.gs[0] = gposei_root*x.gs[0];//new stuff with transformations
	cout<<"x.gs[0]"<<endl<<x.gs[0]<<endl;
  //list of joint angles:
	XmlRpc::XmlRpcValue j_list;
	if(n.getParam("J0", j_list))
	{
		xml2vec(x.r,j_list);
	}
	cout<<"x.r"<<endl<<x.r<<endl;

	if(n.getParam("Jv0", j_list))
	{
		xml2vec(x.dr,j_list);
		cout<<"x.dr"<<endl<<x.dr<<endl;
	}
	mbsmodel->Rec(x, h);

	//Compute Mass matrix:
  MatrixXd M(nb-1 + 6*(!mbsmodel->fixed), nb-1 + 6*(!mbsmodel->fixed));
	mbsmodel->Mass(M,x);
	cout<<"Mass Matrix: "<<endl<<M<<endl;
	//Printing parent indices and joint names and link names:
	for(int count = 0;count <(mbsmodel->nb);count++) 
	{
		cout<<"Links["<<count<<"] "<<mbsmodel->links[count].name<<endl;
		cout<<"pis["<<count<<"]: "<<mbsmodel->pis[count]<<endl;
	}
	Matrix3d rposei_root = gposei_root.topLeftCorner(3,3);
	cout<<"rposei_root"<<endl<<rposei_root<<endl;

	//printing body fixed axes:
	for(int count =0;count<(mbsmodel->nb-1);count++)
	{
		cout<<"Joint body axis["<<mbsmodel->joints[count].name<<"]: "<<mbsmodel->joints[count].a.transpose()<<endl;
		cout<<"x.gs["<<count+1<<"]"<<endl<<x.gs[count+1].topLeftCorner(3,3)<<endl;
		cout<<"mbsmodel->joints["<<count<<"]"<<endl<<mbsmodel->joints[count].gc<<endl;
		cout<<"Joint wrt visual axis["<<(count+1)<<"] "<<(x.gs[count+1].topLeftCorner(3,3)*(mbsmodel->joints[count].gc.topLeftCorner(3,3))*mbsmodel->joints[count].a.head(3)).transpose()<<endl;
	  getchar();
		//cout<<" Dot product: "<<mbsmodel->ag.transpose()*mbsmodel->joints[count].S.head(3)<<endl;

		//cout<<"Joint body axis wrt the visual frame: "<< (rposeroot_i*mbsmodel->joints[count].S.head(3)).transpose()<<endl;
	}

	// initial controls (e.g. hover at one place)
	VectorXd u(mbsmodel->U.n);
	u.setZero();
	if(mbstype == "AIRBASE")
	{
		for(int count = 0;count < nb;count++)
			u[3] += (mbsmodel->links[count].m)*(-mbsmodel->ag(2));
		cout<<"u[3]: "<<u[3]<<endl;
	}
	else if(mbstype == "FLOATBASE")
	{
		for(int count = 0;count < nb;count++)
			u[5] += (mbsmodel->links[count].m)*(-mbsmodel->ag(2));
		cout<<"u[5]: "<<u[5]<<endl;
	}

	//Create states and controls
	xs.resize(N+1,x);
	us.resize(N,u);

	ts.resize(N+1);
	for (int k = 0; k <=N; ++k)
		ts[k] = k*h;

	//Trajectory message initialization
	trajectory.N = N;
	trajectory.statemsg.resize(N+1);
	trajectory.ctrl.resize(N);
	trajectory.time = ts;
	trajectory.rootname = mbsmodel->links[0].name;

	trajectory.statemsg[0].statevector.resize(nb-1);
	trajectory.statemsg[0].names.resize(nb-1);

	for (int i = 0; i < N; ++i) 
	{
		trajectory.statemsg[i+1].statevector.resize(nb-1);
		trajectory.statemsg[i+1].names.resize(nb-1);
		trajectory.ctrl[i].ctrlvec.resize(mbsmodel->U.n);
	}
	//getchar();
	
	// Create timer for iterating	and publishing data
	iteratetimer = n.createTimer(ros::Duration(0.1), simtraj);
	iteratetimer.start();
	//	Dynamic Reconfigure setup Callback ! immediately gets called with default values	
	dynamic_reconfigure::Server<gcop_ctrl::MbsSimInterfaceConfig> server;
	dynamic_reconfigure::Server<gcop_ctrl::MbsSimInterfaceConfig>::CallbackType f;
	f = boost::bind(&paramreqcallback, _1, _2);
	server.setCallback(f);
	ros::spin();
	return 0;
}


	




