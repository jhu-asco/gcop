#ifndef UAV_CAMERA_H
#define UAV_CAMERA_H

#include "opencv/cv.h"
#include "opencv/highgui.h"
#include <iostream>
//#include <pthread.h>
#include "utils.h"
#include <stdint.h>
#include <Eigen/Dense>
#include <fstream>
#define np 5     


namespace gcop {

using namespace Eigen;
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 4, 4> Matrix4d;

  /**
   * Camera Driver
   *
   * Author: Marin Kobilarov -- Copyright (C) 2012
   */
  class Camera {
  public:
    
    Camera(double bl = 1.75, double irows=480, double imcols=640, const double *K=0, const double *D=0);///@bl distance betn two leds; @irows the number of rows in the image @icols the number of columns in the image @K the camera matrix in row major format(9elements) @D the distortion coefficients (plumbob model 5 polynomial coefficients)
    virtual ~Camera();

    bool Pose(Vector6d &q, const cv::Mat &raw, cv::Mat &cont);
		float tolerance;
		int kernel_size;

  protected:
    
    //unsigned int width, height;
    
    cv::Mat filt;    

    cv::Mat rvec;
    cv::Mat tvec;   
    
    double bl;

		//static double _cm[9];//<The camera matrix
		//static double _dc[5];//<Distortion coefficients
		cv::Mat camMatrix;
		cv::Mat dccoeffs;
	  cv::Point2f iprs[np];
		std::ofstream pointfile;
		//bool init_obj;
  };

}

#endif
