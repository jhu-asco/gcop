#ifndef UAV_CAMERA_H
#define UAV_CAMERA_H

#include "opencv/cv.h"
#include "opencv/highgui.h"
#include <iostream>
#include <pthread.h>
#include "utils.h"
#include <stdint.h>
#include <Eigen/Dense>


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
    
    Camera(double bl = 1.75, double irows=480, double imcols=640);
    virtual ~Camera();

    bool Pose(Vector6d &q, const cv::Mat &raw);

  protected:
    
    //unsigned int width, height;
    
    cv::Mat filt;    

    cv::Mat rvec;
    cv::Mat tvec;   
    
    double bl;
  };

}

#endif
