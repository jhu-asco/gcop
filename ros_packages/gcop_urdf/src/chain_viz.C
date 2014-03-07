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
tf::TransformBroadcaster *broadcaster;


// message declarations
geometry_msgs::TransformStamped global_trans;
sensor_msgs::JointState joint_state;
bool recvdmsg = false;

void joint_publish(const gcop_comm::CtrlTraj::ConstPtr& trajectory)
{
	recvdmsg = true;

	//publish zero transform 
	global_trans.header.frame_id = "world";
	global_trans.child_frame_id = "baselink";
	global_trans.header.stamp = ros::Time::now();
	global_trans.transform = trajectory->statemsg[0].basepose;

	broadcaster->sendTransform(global_trans);//send 0 posn


		double pos = 0, dt = 0; //starting wheel position


		for(int i = 1;i<(trajectory->N)+1; i++)
		{
			//joint_state.name.resize(3);
			//joint_state.position.resize(3);
			joint_state.name = trajectory->statemsg[i].names;
			joint_state.header.stamp = ros::Time::now();
			joint_state.position = trajectory->statemsg[i].statevector;
			//send joint state
			joint_pub.publish(joint_state);

			dt = trajectory->time[i] - trajectory->time[i-1];
			ros::Duration(dt).sleep();//sleeps for dt time

			// update transform
			global_trans.header.stamp = ros::Time::now();
			global_trans.transform = trajectory->statemsg[i].basepose;

			//send the transform
			broadcaster->sendTransform(global_trans);
		}
}
int main(int argc, char** argv) {
	ros::init(argc, argv, "state_publisher");
	ros::NodeHandle n;

	broadcaster = new tf::TransformBroadcaster();

	//initializing joint msg
	

	traj_sub = n.subscribe("/mbsdmoc/ctrltraj",1, joint_publish);

	joint_pub = n.advertise<sensor_msgs::JointState>("joint_states", 1);

	ros::Rate loop_rate(100);
	while(ros::ok())
	{
		if(recvdmsg)
		{
			//publish zero transform 
			joint_state.header.stamp = ros::Time::now();
			global_trans.header.stamp = ros::Time::now();
			joint_pub.publish(joint_state);
			broadcaster->sendTransform(global_trans);
		}
		ros::spinOnce();
		loop_rate.sleep();
	}


	
		return 0;
}
