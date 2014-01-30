/** This node takes the ctrl trajectory data for visualization. It publishes the joint message
* to the robot state publisher and also publishes the baselink transformation needed to move the robot.
*/
#include "ros/ros.h"
#include <sensor_msgs/JointState.h>
#include <tf/transform_broadcaster.h>
#include <visualization_msgs/Marker.h>
#include "gcop_comm/CtrlTraj.h"


ros::Subscriber traj_sub;
ros::Publisher joint_pub;
ros::Publisher finaljoint_pub;
ros::Publisher traj_pub;
ros::Publisher goal_pub;
tf::TransformBroadcaster *broadcaster;


// message declarations
geometry_msgs::TransformStamped global_trans;
geometry_msgs::TransformStamped finalgoal_trans;
sensor_msgs::JointState joint_state;
sensor_msgs::JointState finaljoint_state;
visualization_msgs::Marker line_strip;
visualization_msgs::Marker goal;

void joint_publish(const gcop_comm::CtrlTraj::ConstPtr& trajectory)
{
	//This function is specific to car and should be changed for other vehicles
	//ros param define stuff	
	bool animate = false;
	ros::param::getCached("/animate", animate);

	goal.header.frame_id = "world";
	goal.header.stamp = ros::Time::now();
	goal.ns ="goalmarker";
	goal.id = 2;
	goal.type = visualization_msgs::Marker::ARROW;
	goal.action = visualization_msgs::Marker::ADD;
	goal.pose.position.x = trajectory->finalgoal.statevector[0];
	goal.pose.position.y = trajectory->finalgoal.statevector[1];
	goal.pose.position.z = 0.1;
	goal.pose.orientation = tf::createQuaternionMsgFromYaw(trajectory->finalgoal.statevector[2]);
	goal.scale.x = 1;
	goal.scale.y = 1;
	goal.scale.z = trajectory->finalgoal.statevector[3];

	goal.color.r = 0.8;
	goal.color.g = 0.0;
	goal.color.b = 0.0;
	goal.color.a = 1.0;

	goal_pub.publish(goal);

	//publish zero transform 
	global_trans.header.stamp = ros::Time::now();
	global_trans.transform.translation.x = trajectory->statemsg[0].statevector[0];
	global_trans.transform.translation.y = trajectory->statemsg[0].statevector[1];
	global_trans.transform.translation.z = .1;
	global_trans.transform.rotation = tf::createQuaternionMsgFromYaw(trajectory->statemsg[0].statevector[2]);

	broadcaster->sendTransform(global_trans);//send 0 posn

	finalgoal_trans.header.stamp = ros::Time::now();
	finalgoal_trans.transform.translation.x = trajectory->finalgoal.statevector[0];
	finalgoal_trans.transform.translation.y = trajectory->finalgoal.statevector[1];
	finalgoal_trans.transform.translation.z = .1;
	finalgoal_trans.transform.rotation = tf::createQuaternionMsgFromYaw(trajectory->finalgoal.statevector[2]);

	broadcaster->sendTransform(finalgoal_trans);//send 0 posn

	finaljoint_state.header.stamp = ros::Time::now();
	finaljoint_state.position[2] = 0;
	//steering angle
	finaljoint_state.position[0] = trajectory->ctrl[trajectory->N-1].ctrlvec[1];
	finaljoint_state.position[1] = trajectory->ctrl[trajectory->N-1].ctrlvec[1];
	//send joint state
	finaljoint_pub.publish(finaljoint_state);



	//just publish the trajectory
	line_strip.header.stamp  = ros::Time::now();
	line_strip.points.resize(trajectory->N + 1);

	for(int i =0;i<(trajectory->N)+1; i++)
	{
		//geometry_msgs::Point p;
		line_strip.points[i].x = trajectory->statemsg[i].statevector[0];
		line_strip.points[i].y =  trajectory->statemsg[i].statevector[1];
		line_strip.points[i].z = 0.1;
		//line_strip.points.push_back(p);
	}
	traj_pub.publish(line_strip);


	if(animate)
	{
		double pos = 0, dt = 0; //starting wheel position


		for(int i = 1;i<(trajectory->N)+1; i++)
		{
			// convert prev wheel vel to new pos
			//joint_state.position[0] = pos + trajectory->ctrl[i].ctrlvec[0]*dt ;
			//for now pos[0] = 0;
			joint_state.header.stamp = ros::Time::now();
			joint_state.position[2] = 0;
			//steering angle
			joint_state.position[0] = trajectory->ctrl[i-1].ctrlvec[1];
			joint_state.position[1] = trajectory->ctrl[i-1].ctrlvec[1];
			//send joint state
			joint_pub.publish(joint_state);

			dt = trajectory->time[i] - trajectory->time[i-1];
			ros::Duration(dt).sleep();//sleeps for dt time

			// update transform
			global_trans.header.stamp = ros::Time::now();
			global_trans.transform.translation.x = trajectory->statemsg[i].statevector[0];
			global_trans.transform.translation.y = trajectory->statemsg[i].statevector[1];
			global_trans.transform.translation.z = .1;
			global_trans.transform.rotation = tf::createQuaternionMsgFromYaw(trajectory->statemsg[i].statevector[2]);

			//send the transform
			broadcaster->sendTransform(global_trans);
		}
	}
}
int main(int argc, char** argv) {
	ros::init(argc, argv, "state_publisher");
	ros::NodeHandle n;

	broadcaster = new tf::TransformBroadcaster();

	//initializing joint msg
	joint_state.name.resize(3);
	joint_state.header.frame_id = "movingcar";
	joint_state.position.resize(3);
	joint_state.name[0] = "base_to_frontwheel1";
	joint_state.name[1] = "base_to_frontwheel2";
	joint_state.name[2] = "base_to_backwheel1";

	//initializing finaljoint msg
	finaljoint_state.name.resize(3);
	finaljoint_state.header.frame_id = "goalcar";
	finaljoint_state.position.resize(3);
	finaljoint_state.name[0] = "base_to_frontwheel1";
	finaljoint_state.name[1] = "base_to_frontwheel2";
	finaljoint_state.name[2] = "base_to_backwheel1";


	traj_sub = n.subscribe("/dmoc/ctrltraj",1, joint_publish);

	joint_pub = n.advertise<sensor_msgs::JointState>("/movingcar/joint_states", 1);
	finaljoint_pub = n.advertise<sensor_msgs::JointState>("/goalcar/joint_states", 1);
	traj_pub = n.advertise<visualization_msgs::Marker>("desired_traj", 1);
	goal_pub = n.advertise<visualization_msgs::Marker>("goal_config", 1);
	ros::Publisher marker_pub = n.advertise<visualization_msgs::Marker>("Obstacle_marker", 1);

	ros::Rate loop_rate(100);
	//initialize the position and state
	global_trans.header.frame_id = "world";
	global_trans.child_frame_id = "/movingcar/baselink";
	global_trans.transform.translation.x = 0;
	global_trans.transform.translation.y = 0;
	global_trans.transform.translation.z = .1;
	global_trans.transform.rotation = tf::createQuaternionMsgFromYaw(0);

	finalgoal_trans.header.frame_id = "world";
	finalgoal_trans.child_frame_id = "/goalcar/baselink";
	finalgoal_trans.transform.translation.x = 0;
	finalgoal_trans.transform.translation.y = 0;
	finalgoal_trans.transform.translation.z = .1;
	finalgoal_trans.transform.rotation = tf::createQuaternionMsgFromYaw(0);


	joint_state.position[0] = 0;
	joint_state.position[1] = 0;
	joint_state.position[2] = 0; //default positions

	finaljoint_state.position[0] = 0;
	finaljoint_state.position[1] = 0;
	finaljoint_state.position[2] = 0; //default positions


	line_strip.header.frame_id = "/world";
	line_strip.ns = "traj";
	line_strip.action = visualization_msgs::Marker::ADD;
	line_strip.pose.orientation.w = 1.0;
	line_strip.id = 1;
	line_strip.type = visualization_msgs::Marker::LINE_STRIP;
	line_strip.scale.x = 0.1;
	line_strip.color.b = 1.0;
	line_strip.color.a = 1.0;

	//Adding HardCoded Obstacle:
	visualization_msgs::Marker marker;
	marker.ns = "cylinder_obst";
	marker.id = 0;
	marker.type = visualization_msgs::Marker::CYLINDER;
	marker.action = visualization_msgs::Marker::ADD;

	marker.header.frame_id = "/world";
	marker.pose.position.x = -0.5;
	marker.pose.position.y = -0.5;
	marker.pose.position.z = 0.0;
	marker.pose.orientation = tf::createQuaternionMsgFromYaw(0);

	// Set the scale of the marker -- 1x1x1 here means 1m on a side
	marker.scale.x = 0.1;
	marker.scale.y = 0.1;
	marker.scale.z = 0.25;

	// Set the color -- be sure to set alpha to something non-zero!
	marker.color.r = 0.0f;
	marker.color.g = 1.0f;
	marker.color.b = 0.0f;
	marker.color.a = 1.0;

	marker.lifetime = ros::Duration();

	while(ros::ok())
	{
		marker_pub.publish(marker);
		//publish the last trankfor.  and joint position
		joint_state.header.stamp = ros::Time::now();
		joint_pub.publish(joint_state);
		//finaljoints
		finaljoint_state.header.stamp = ros::Time::now();
		finaljoint_pub.publish(finaljoint_state);
		//transform to movingcar
		global_trans.header.stamp = ros::Time::now();
		broadcaster->sendTransform(global_trans);//send 0 posn
		//transform to finalgoal 
		finalgoal_trans.header.stamp = ros::Time::now();
		broadcaster->sendTransform(finalgoal_trans);//send finalgoal 
		//Marker
		marker.header.stamp = ros::Time::now();
		ros::spinOnce();
		loop_rate.sleep();
	}

	return 0;
}
