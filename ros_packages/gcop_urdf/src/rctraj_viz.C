/** This node takes the ctrl trajectory data for visualization. It publishes the joint message
* to the robot state publisher and also publishes the baselink transformation needed to move the robot.
*/
#include "ros/ros.h"
#include <sensor_msgs/JointState.h>
#include <tf/transform_broadcaster.h>
#include "gcop_comm/CtrlTraj.h"


ros::Subscriber traj_sub;
ros::Publisher joint_pub;
tf::TransformBroadcaster *broadcaster;

// message declarations
geometry_msgs::TransformStamped global_trans;
sensor_msgs::JointState joint_state;

void joint_publish(const gcop_comm::CtrlTraj::ConstPtr& trajectory)
{
	//This function is specific to car and should be changed for other vehicles
	
	double pos = 0, dt = 0; //starting wheel position

	//publish zero transform 
		global_trans.header.frame_id = "world";
		global_trans.child_frame_id = "baselink";
		global_trans.header.stamp = ros::Time::now();
		global_trans.transform.translation.x = trajectory->statemsg[0].statevector[0];
		global_trans.transform.translation.y = trajectory->statemsg[0].statevector[1];
		global_trans.transform.translation.z = .1;
		global_trans.transform.rotation = tf::createQuaternionMsgFromYaw(trajectory->statemsg[0].statevector[2]);

		broadcaster->sendTransform(global_trans);//send 0 posn


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
int main(int argc, char** argv) {
	ros::init(argc, argv, "state_publisher");
	ros::NodeHandle n;

	broadcaster = new tf::TransformBroadcaster();

	//initializing joint msg
	joint_state.name.resize(3);
	joint_state.position.resize(3);
	joint_state.name[0] = "base_to_frontwheel1";
	joint_state.name[1] = "base_to_frontwheel2";
	joint_state.name[2] = "base_to_backwheel1";


	traj_sub = n.subscribe("/dmoc/ctrltraj",1, joint_publish);

	joint_pub = n.advertise<sensor_msgs::JointState>("joint_states", 1);

	ros::Rate loop_rate(100);
	//initialize the position and state
	global_trans.header.frame_id = "world";
	global_trans.child_frame_id = "baselink";
	global_trans.transform.translation.x = 0;
	global_trans.transform.translation.y = 0;
	global_trans.transform.translation.z = .1;
	global_trans.transform.rotation = tf::createQuaternionMsgFromYaw(0);

	joint_state.position[0] = 0;
	joint_state.position[1] = 0;
	joint_state.position[2] = 0; //default positions


	while(ros::ok())
	{
		//publish the last transform  and joint position
		joint_state.header.stamp = ros::Time::now();
		joint_pub.publish(joint_state);
		global_trans.header.stamp = ros::Time::now();
		broadcaster->sendTransform(global_trans);//send 0 posn
		ros::spinOnce();
		loop_rate.sleep();
	}

	return 0;
}
