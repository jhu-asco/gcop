#include "ros/ros.h"
#include <iomanip>
#include <iostream>
#include <dynamic_reconfigure/server.h>
#include "gcop_comm/CtrlTraj.h"//msg for publishing ctrl trajectory
#include <gcop/urdf_parser.h>
#include "tf/transform_datatypes.h"
#include <gcop/se3.h>

using namespace std;
using namespace Eigen;
using namespace gcop;

//ros messages
gcop_comm::CtrlTraj trajectory;

//Publisher
ros::Publisher trajpub;

void q2transform(geometry_msgs::Transform &transformmsg, Vector6d &bpose)
{
	tf::Quaternion q;
	q.setEulerZYX(bpose[2],bpose[1],bpose[0]);
	tf::Vector3 v(bpose[3],bpose[4],bpose[5]);
	tf::Transform tftransform(q,v);
	tf::transformTFToMsg(tftransform,transformmsg);
	//cout<<"---------"<<endl<<transformmsg.position.x<<endl<<transformmsg.position.y<<endl<<transformmsg.position.z<<endl<<endl<<"----------"<<endl;
}

int main(int argc, char** argv)
{
	ros::init(argc, argv, "chainload");
	ros::NodeHandle n;
	//Initialize publisher
	trajpub = n.advertise<gcop_comm::CtrlTraj>("ctrltraj",2);
	//get parameter for xml_string:
	string xml_string, xml_filename;
	if(!ros::param::get("/robot_description", xml_string))
	{
		ROS_ERROR("Could not fetch xml file name");
		return 0;
	}
	//cout<<xml_filename<<endl;
	//fstream xml_file(xml_filename.c_str(), fstream::in);
	/*while ( xml_file.good() )
	{
		string line;
		getline( xml_file, line);
		xml_string += (line + "\n");
	}
	xml_file.close();
	*/

	string mbstype;
	n.getParam("basetype",mbstype);
	Matrix4d gposei_root;
	boost::shared_ptr<Mbs> mbsmodel = gcop_urdf::mbsgenerator(xml_string,gposei_root,mbstype);
	//set no gravity:
	mbsmodel->ag << 0, 0, -9.81;

	//check if mbsmodel is good by using it.
	//first try printing the mbsmodel params:
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
	//Looks ok till now
	//Using it:

  double h = .01;
  int N = 1000;     // discrete trajectory segments
  double tf = h*N;   // time-horizon
	vector<double> ts;
	ts.resize(N+1);
	int nb = mbsmodel->nb;

	for (int k = 0; k <=N; ++k)
		ts[k] = k*h;
  //Trajectory message initialization
	trajectory.N = N;
	trajectory.statemsg.resize(N+1);
	trajectory.ctrl.resize(N);
	trajectory.time = ts;
	//trajectory.finalgoal.statevector.resize(4);

  //mbsmodel->debug = true;


  MbsState x(nb);
  x.gs[0].setIdentity();
	x.gs[0](0,3) = 1.0;
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

	//cout<<"M_PI"<<endl;
	//cout<<M_PI<<endl;

  // initial controls (e.g. hover at one place)
  VectorXd u(6+ (nb-1));
  u.setZero();
	u[5] = 30*9.81;
  vector<VectorXd> us(N, u);

  //struct timeval timer;
  //  dmoc.debug = false; // turn off debug for speed

	//Initialize root pose
	//Initialize 0 state
	//cout<<"----"<<endl<<xs[0].gs[0]<<"----"<<endl;

	while(ros::ok())
	{
		Vector6d bpose;
		gcop::SE3::Instance().g2q(bpose, xs[0].gs[0]);
		q2transform(trajectory.statemsg[0].basepose,bpose);
		trajectory.statemsg[0].statevector.resize(nb-1);
		trajectory.statemsg[0].names.resize(nb-1);
		for(int count1 = 0;count1 < nb-1;count1++)
		{
			trajectory.statemsg[0].statevector[count1] = xs[0].r[count1];
			trajectory.statemsg[0].names[count1] = mbsmodel->joints[count1].name;
		}
		//cout<<"xs0"<<xs[0].r[0]<<"\t"<<xs[0].r[1]<<"\t"<<trajectory.statemsg[0].statevector[0]<<trajectory.statemsg[0].statevector[1]<<"\t"<<endl;

		for (int i = 0; i < N; ++i) 
		{
			//timer_start(timer);
			mbsmodel->Step(xs[i+1], i*h, xs[i], us[i], h);
			//long te = timer_us(timer);
			//cout << "Iteration #" << i << ": took " << te << " us." << endl;
			trajectory.statemsg[i+1].statevector.resize(nb-1);
			trajectory.statemsg[i+1].names.resize(nb-1);
			gcop::SE3::Instance().g2q(bpose, xs[i+1].gs[0]);
			q2transform(trajectory.statemsg[i+1].basepose,bpose);
			for(int count1 = 0;count1 < nb-1;count1++)
			{
				trajectory.statemsg[i+1].statevector[count1] = xs[i+1].r[count1];
				trajectory.statemsg[i+1].names[count1] = mbsmodel->joints[count1].name;
			}
			trajectory.ctrl[i].ctrlvec.resize(6+nb-1);
			for(int count1 = 0;count1 < 6+nb-1;count1++)
			{
				trajectory.ctrl[i].ctrlvec[count1] = us[i](count1);
			}
		}
		//final goal:
		/*for(int count = 0;count<4;count++)
			{
			trajectory.finalgoal.statevector[count] = xf(count);
			}
		 */
		//cout<<"xs0"<<xs[0].r[0]<<"\t"<<xs[0].r[1]<<"\t"<<trajectory.statemsg[0].statevector[0]<<trajectory.statemsg[0].statevector[1]<<"\t"<<endl;
		trajpub.publish(trajectory);
		ros::spinOnce();

		cout << "dr=" << xs[1].dr << endl;

		cout << "done!" << endl;
		ros::Duration(0.05).sleep();
		getchar();
	}

	return 0;
}


	//boost::shared_ptr<gcop_urdf::ModelInterface> robot = gcop_urdf::parseURDF(xml_string);
	/*
	if (!robot){
		cerr << "ERROR: Model Parsing the xml failed" << endl;
		return -1;
	}
	cout << "robot name is: " << robot->getName() << endl;

	// get info from parser
	cout << "---------- Successfully Parsed XML ---------------" << endl;
	// get root link
	boost::shared_ptr<const Link> root_link=robot->getRoot();
	if (!root_link) return -1;

	cout << "root Link: " << root_link->name << " has " << root_link->child_links.size() << " child(ren)" << endl;
	*/




