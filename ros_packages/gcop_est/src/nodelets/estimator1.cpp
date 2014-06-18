#include <ros/ros.h>
#include <nodelet/nodelet.h>
#include <image_transport/image_transport.h>
#include <cv_bridge/cv_bridge.h>
#include <Eigen/Dense>
#include <gcop/so3.h>
#include <gcop/camera.h>
#include <iostream>
#include <tf/transform_broadcaster.h>

/* This nodelet is a wrapper for camera library (name is a misnomer) in gcop. It finds the 6dof pose of the camera with respect to an known object.
	 The code has mostly been copied from Rectify.cpp of image_proc ros package
 */

namespace gcop_est {
	using namespace Eigen;
	using namespace gcop;
	using namespace std;
	class PoseEstimator1 : public nodelet::Nodelet
	{
		// ROS communication
		boost::shared_ptr<image_transport::ImageTransport> it_;
		image_transport::CameraSubscriber sub_camera_;

		//boost::mutex connect_mutex_;

		//ros::Publisher pose_publisher;
		bool initcamera;
		int length_object;
		boost::shared_ptr<gcop::Camera> cam;
		Vector6d q;
		Matrix3d g;
		//Vector3d rpy;//We will change this as 	quaternion for later use
		Vector4d wxyz;
		//Vector3d xyz;

		//Member Functions:
		virtual void onInit();
		//void connectCb();
		void imageCb(const sensor_msgs::ImageConstPtr& image_msg,
				const sensor_msgs::CameraInfoConstPtr& info_msg);
	};//Class PoseEstimator1
	void PoseEstimator1::onInit()
	{
		ros::NodeHandle &nh         = getNodeHandle();
		ros::NodeHandle &private_nh = getPrivateNodeHandle();
		initcamera = false;//gcop camera class is uninitialized
		length_object = 3.5;//Hardcoded for now
		it_.reset(new image_transport::ImageTransport(nh));
		//Declare ros publisher for pose or use tf
		// Monitor whether anyone is subscribed to the output
		//image_transport::SubscriberStatusCallback connect_cb = boost::bind(&RectifyNodelet::connectCb, this);
		image_transport::TransportHints hints("raw", ros::TransportHints(), private_nh);
		sub_camera_ = it_->subscribeCamera("image", 1, &PoseEstimator1::imageCb, this, hints);
	}

	// Handles (un)subscribing when clients (un)subscribe
	/*void RectifyNodelet::connectCb()
	{
		boost::lock_guard<boost::mutex> lock(connect_mutex_);
		if (pub_rect_.getNumSubscribers() == 0)
			sub_camera_.shutdown();
		else if (!sub_camera_)
		{
			image_transport::TransportHints hints("raw", ros::TransportHints(), getPrivateNodeHandle());
			sub_camera_ = it_->subscribeCamera("image_mono", queue_size_, &RectifyNodelet::imageCb, this, hints);
		}
	}*/
	//Main function where image is received
	void PoseEstimator1::imageCb(const sensor_msgs::ImageConstPtr& image_msg,
			const sensor_msgs::CameraInfoConstPtr& info_msg)
	{
		static tf::TransformBroadcaster br;
		tf::Transform result_transform;
		// Verify camera is actually calibrated
		if (info_msg->K[0] == 0.0) {
			cout<<"Camera not calibrated"<<endl;
			return;
		}
		if(!initcamera)
		{
			cam.reset(new gcop::Camera(length_object,info_msg->height,info_msg->width));
			initcamera = true;//Camera class is initialized
		}
		// Create cv::Mat views onto both buffers
		const cv::Mat image = cv_bridge::toCvShare(image_msg)->image;
		if(cam)//If cam is initialized
		{
			//Use the camera library to process for pose
			bool ok = cam->Pose(q,image);//Find the pose
			if (!ok) {
				ROS_WARN("No Object");
				return;
			}
			//xyz =q.tail<3>();//First three are posn
			SO3::Instance().exp(g,q.head<3>());//Convert the exponential coords to matrix
			SO3::Instance().g2quat(wxyz,g);//Convert g to rpy
			result_transform.setOrigin(tf::Vector3(q(3),q(4),q(5)));
			result_transform.setRotation(tf::Quaternion(wxyz(1),wxyz(2),wxyz(3),wxyz(0)));
			//br.sendTransform(tf::StampedTransform(result_transform,image_msg->header.stamp,"camera","object"));//transform from camera(parent) frame to object(child) frame
			br.sendTransform(tf::StampedTransform(result_transform,ros::Time::now(),"camera","object"));//transform from camera(parent) frame to object(child) frame
			/*
			cout<<"Q output: "<<q.transpose()<<endl;

			cout << "rpy=" << rpy.transpose() << "xyz=" << xyz.transpose() << endl;

			cout << "x=" << -xyz[:2] << " y=" << xyz[0] << " yaw=" << -rpy[0]*(180.0/M_PI) << endl;//Will move this into tf publishing TODO
			*/
		}
	}
}//namespace gcop_est

	// Register nodelet
#include <pluginlib/class_list_macros.h>
	PLUGINLIB_DECLARE_CLASS(gcop_est, estimator1, gcop_est::PoseEstimator1, nodelet::Nodelet)

