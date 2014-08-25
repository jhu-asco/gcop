#include "camera.h"
#include "utils.h"
#include <iostream>
#include "so3.h"
#include <dc1394/dc1394.h>
#include <signal.h>

using namespace std;
using namespace gcop;
using namespace Eigen;
//using namespace itpp;

//Separating image acquisition from the library 
// This way we can even work on offline images
    dc1394camera_t *camera = 0;
    dc1394video_frame_t *frame = NULL;
    //dc1394featureset_t features = 0;
    dc1394_t * d = 0;
    dc1394camera_list_t * list = 0;
    dc1394error_t err;

    cv::Mat raw; //Raw Image being set from the camera
    cv::Mat cont; //Raw Image being set from the camera


void signal_callback_handler(int signum)//Terminate signal
{
	cout<<"Stopping Camera ........."<<endl;
  dc1394_video_set_transmission(camera, DC1394_OFF);
  dc1394_capture_stop(camera);
  dc1394_camera_free(camera);
	exit(0);
}

void cleanup_and_exit(dc1394camera_t *camera)
{
    dc1394_video_set_transmission(camera, DC1394_OFF);
    dc1394_capture_stop(camera);
    dc1394_camera_free(camera);
    exit(1);
}

int StartCamera()
{
  if (camera) 
		signal_callback_handler(0);//Stop Camera if it is already started
  
  d = dc1394_new ();
  if (!d)
    return 1;
  err=dc1394_camera_enumerate (d, &list);
  DC1394_ERR_RTN(err,"Failed to enumerate cameras");
  
  if (list->num == 0) {
    dc1394_log_error("No cameras found");
    return 1;
  }
  
  camera = dc1394_camera_new (d, list->ids[0].guid);
  if (!camera) {
    dc1394_log_error("Failed to initialize camera with guid %llx", list->ids[0].guid);
    return 1;
  }
  dc1394_camera_free_list (list);
  
  //    printf("Using camera with GUID %l\n", camera->guid);
  
  
  /*-----------------------------------------------------------------------
   *  setup capture
   *-----------------------------------------------------------------------*/
  
  //    err=dc1394_video_set_iso_speed(camera, DC1394_ISO_SPEED_400);
  //    DC1394_ERR_CLN_RTN(err,cleanup_and_exit(camera),"Could not set iso speed");
  
  //err=dc1394_video_set_mode(camera, DC1394_VIDEO_MODE_FORMAT7_0);
  //    err=dc1394_video_set_mode(camera, DC1394_VIDEO_MODE_1280x960_MONO8);
  //    
  err=dc1394_video_set_mode(camera, DC1394_VIDEO_MODE_640x480_MONO8);
  DC1394_ERR_CLN_RTN(err,cleanup_and_exit(camera),"Could not set video mode\n");
  
  
  //    err=dc1394_video_set_framerate(camera, DC1394_FRAMERATE_15);
  //    DC1394_ERR_CLN_RTN(err,cleanup_and_exit(camera),"Could not set framerate\n");
  
    unsigned int bmin, bmax;
    
     dc1394_feature_set_mode(camera, DC1394_FEATURE_GAIN, DC1394_FEATURE_MODE_MANUAL);
     dc1394_feature_get_boundaries(camera,DC1394_FEATURE_GAIN, &bmin, &bmax);     
     dc1394_feature_set_value(camera,DC1394_FEATURE_GAIN, bmin);
     //    dc1394_feature_set_value(camera,DC1394_FEATURE_GAIN,160);
    

    dc1394_feature_set_mode(camera, DC1394_FEATURE_FRAME_RATE, DC1394_FEATURE_MODE_MANUAL);
    dc1394_feature_set_value(camera,DC1394_FEATURE_FRAME_RATE,30);
    
    
    //    dc1394_feature_set_value(camera,DC1394_FEATURE_SHUTTER,580);
    dc1394_feature_set_mode(camera, DC1394_FEATURE_SHUTTER, DC1394_FEATURE_MODE_MANUAL);
     dc1394_feature_get_boundaries(camera,DC1394_FEATURE_SHUTTER, &bmin, &bmax);     
     //    dc1394_feature_set_value(camera,DC1394_FEATURE_SHUTTER, bmax);
    dc1394_feature_set_value(camera,DC1394_FEATURE_SHUTTER, bmin);
    

    dc1394_feature_set_value(camera,DC1394_FEATURE_BRIGHTNESS,1);
    dc1394_feature_set_mode(camera, DC1394_FEATURE_BRIGHTNESS, DC1394_FEATURE_MODE_MANUAL);

    dc1394_feature_set_value(camera,DC1394_FEATURE_EXPOSURE,1);
    dc1394_feature_set_mode(camera, DC1394_FEATURE_EXPOSURE, DC1394_FEATURE_MODE_MANUAL);

    dc1394_feature_set_value(camera,DC1394_FEATURE_GAMMA,512);
    dc1394_feature_set_mode(camera, DC1394_FEATURE_GAMMA, DC1394_FEATURE_MODE_MANUAL);

    //err=dc1394_capture_setup(camera, 10, DC1394_CAPTURE_FLAGS_DEFAULT);
    err=dc1394_capture_setup(camera,1, 0);
    DC1394_ERR_CLN_RTN(err,cleanup_and_exit(camera),"Could not setup camera-\nmake sure that the video mode and framerate are\nsupported by your camera\n");

       

    /*-----------------------------------------------------------------------
     *  have the camera start sending us data
     *-----------------------------------------------------------------------*/
    err=dc1394_video_set_transmission(camera, DC1394_ON);
    DC1394_ERR_CLN_RTN(err,cleanup_and_exit(camera),"Could not start camera iso transmission\n");

    /*-----------------------------------------------------------------------
     *  capture one frame
     *-----------------------------------------------------------------------*/
    err = dc1394_capture_dequeue(camera, DC1394_CAPTURE_POLICY_WAIT, &frame);
    DC1394_ERR_CLN_RTN(err,cleanup_and_exit(camera),"Could not capture a frame\n");
   
    /*
    cap.set(FlyCapture2::BRIGHTNESS,0);
    cap.set(FlyCapture2::AUTO_EXPOSURE, 1);

    //    cap.set(FlyCapture2::WHITE_BALANCE, 1);
    cap.set(FlyCapture2::GAMMA, 512);
    cap.set(FlyCapture2::SHUTTER, 1077);
    cap.set(FlyCapture2::GAIN, 160);
    //    cap.set(FlyCapture2::SHARPNESS, 160);
    */
    raw = cv::Mat(frame->size[1], frame->size[0], CV_8UC1, frame->image);
		return 0;
}

bool CaptureCamera()
{
  dc1394_capture_enqueue(camera, frame);
  
  err = dc1394_capture_dequeue(camera, DC1394_CAPTURE_POLICY_WAIT, &frame);
  DC1394_ERR_CLN_RTN(err,cleanup_and_exit(camera),"Could not capture a frame\n");
    
  memcpy(raw.data, frame->image, frame->size[0]*frame->size[1]*sizeof(char));
  return true;
}

int main(int argc, char** argv)
{
	// Register signal and signal handler
	signal(SIGINT, signal_callback_handler);
	//Initialize Camera
	StartCamera();
	Camera cam(3.5, frame->size[1], frame->size[0]);//Just as given we will also print them to see
  cout<<"Frame size: "<<frame->size[1]<<"\t"<<frame->size[0]<<endl;

	//  mat gp = eye(3);
	//  SE2M::Instance().exp(gp, vec( (double[3]){M_PI,0,0}, 3));
	//  mat g = eye(3);

	Vector6d q;
	Matrix3d g;
	Vector3d rpy;
	Vector3d xyz;
	SO3 so3 = SO3::Instance();



	//mat R(3,3);
	//vec x(3);

	while (1) {
		CaptureCamera();//Get a raw frame  from camera
		bool ok = cam.Pose(q,raw,cont);//Compute the pose
		if (!ok) {
			cout << "no object" << endl;
			continue;
		} 

		//    cout << "q= " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " " << q[4] << " " << q[5] << endl;
		//      cout << "yaw= " << q[2] << " x=" << q[5] << " " << q[3] << endl;

		//R = SO3::exp(vec(q, 3));
		xyz =q.tail<3>();//First three are posn
		so3.exp(g,q.head<3>());//Convert the exponential coords to matrix
		so3.g2q(rpy,g);//Convert g to rpy
		cout<<"Q output: "<<q.transpose()<<endl;

		cout << "rpy=" << rpy.transpose() << "xyz=" << xyz.transpose() << endl;

		cout << "x=" << -xyz[2] << " y=" << xyz[0] << " yaw=" << -rpy[0]*(180.0/M_PI) << endl;
	}
}
