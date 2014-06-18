#include "camera.h"
#include <stdio.h>
#include <string.h>
#include <limits>

/* This library works on 1394 firewire cameras to find the 6DOF pose
 * of the camera with respect to a known engineered object
 * The camera needs to be calibrated and also the object lengths
 * should be know before hand 
 */

using namespace cv;
using namespace gcop;
using namespace Eigen;
using namespace std;

static double area(const Point &a, const Point &b, const Point &c)
{
  Point v = b - a;
  Point w = c - a;
  return v.x*w.y - v.y*w.x;  
}

static bool inters(const Point &a, const Point &b, const Point &c, const Point &d)
{
  return (area(a, b, c)*area(a, b, d) < 0) && (area(c, d, a)*area(c, d, b) < 0);
}

static void swap(Point &a, Point &b)
{
  Point t = a;
  a = b;
  b = t;
}


/*
E = B-A = ( Bx-Ax, By-Ay )
F = D-C = ( Dx-Cx, Dy-Cy ) 
P = ( -Ey, Ex )
h = ( (A-C) * P ) / ( F * P )
*/


Camera::Camera(double bl, double imrows, double imcols)
{
	this->bl = bl;
	filt = cv::Mat(imrows, imcols, CV_8UC1);
	rvec = cv::Mat::zeros(3,1, CV_64FC1);
	tvec = cv::Mat::zeros(3,1, CV_64FC1);
	tvec.at<float>(2,0) = 10;
}

Camera::~Camera()
{
}



bool Camera::Pose(Vector6d &q, const cv::Mat &raw)//Input Raw cv image and output 6dof pose
{
  threshold(raw, filt, 150, 255, THRESH_BINARY);
  blur( filt, filt, Size(3,3) );
  
  vector< vector<Point> > contours;
  findContours(filt, contours, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);
  Scalar color( 50, 50, 50 );
  drawContours(raw, contours, -1, color,CV_FILLED);

  const int np = 5;     
  /*
    b a c
    
    d   e
  */
  
#define FLIP_AXIS
  
  //      float r = 4.5;
  //  float bl = 3.5;
#ifdef FLIP_AXIS
  
  
  static float _oprs[np][3] = {{(float)bl, 0, 0},//Need to provide a way to access these from outside the class TODO 
                        {(float)bl, -(float)bl, 0},
                        {(float)bl, (float)bl, 0},
                        {-(float)bl, -(float)bl, 0}, 
                        {-(float)bl, (float)bl, 0}};
#else
  
  static float _oprs[np][3] = {{0, (float)bl, 0},
                               {-(float)bl, (float)bl, 0},
                               {(float)bl, (float)bl, 0},
                               {-(float)bl, -(float)bl, 0}, 
                               {(float)bl, -(float)bl, 0}};
#endif
  
  Point2f iprs[np];
  
  Mat O(np, 3, CV_32FC1,_oprs);
  Mat I(np, 2, CV_32FC1);
  
  vector< vector<Point> >::iterator ci;
  if (contours.size() >= np) {
    int i = 0;      
    for (ci = contours.begin(); ci != contours.end(); ++ci) { 
      vector<Point>::iterator pi;
      double x = 0, y = 0;
      for (pi = ci->begin(); pi != ci->end(); ++pi) {
        x += pi->x;
        y += pi->y;
      }
      x /= (*ci).size();
      y /= (*ci).size();
      
      if (i < np) {
        iprs[i].x = x;
        iprs[i].y = y;
        ++i;
      }
    }
    
    double min = numeric_limits<double>::max();

    int cs[np];
    for (int i = 0; i < np; ++i) {
      for (int j = 0; j < np; ++j) {
        if (j == i)
          continue;
        for (int k = 0; k < np; ++k) {
          if (k == i || k == j)
            continue;
          
          Point2f& a = iprs[i];
          Point2f& b = iprs[j];        
          Point2f& c = iprs[k];
          
          Point2f ba = b - a;
          Point2f ca = c - a;
          double ban = sqrt(ba.x*ba.x + ba.y*ba.y);
          double can = sqrt(ca.x*ca.x + ca.y*ca.y);
          // add slopes -- middle point should have minimum value
          double r = ba.dot(ca)/(ban*can);
          if (r < min) {
            min = r;
            cs[0] = i;
            cs[1] = j;
            cs[2] = k;
          }              
        }
      }
    }
    
    // mask the first three
    bool ms[5] = {0,0,0,0,0};
    ms[cs[0]] = 1;
    ms[cs[1]] = 1;
    ms[cs[2]] = 1;
    
    // fill in the remaining two
    for (int i = 0, j = 3; i < 5; ++i) {
      if (!ms[i]) {
        cs[j] = i;
        ++j;
      }
    }
    
    Point2f& a = iprs[cs[0]];
    Point2f& b = iprs[cs[1]];        
    Point2f& c = iprs[cs[2]];
    Point2f& d = iprs[cs[3]];
    Point2f& e = iprs[cs[4]];
        
    if (inters(b, d, c, e)) {
      swap(d, e);
    }
    
    if (area(a,b,d) > 0) {
      swap(b,c);
      swap(d,e);
    }
    
    //        CvFont font = InitFont(CV_FONT_HERSHEY_SIMPLEX, 1, 1);
    
    for (int i = 0; i < np; ++i) {
      I.at<float>(i,0) = iprs[cs[i]].x;
      I.at<float>(i,1) = iprs[cs[i]].y;
    }
    
    // firefly MTC
    static double _cm[9] = {623.77133283076, 0, 325.476798682001,
														 0, 624.840987895719, 180.290150726724,
														  0, 0, 1};
    /* scorpion
    static double _cm[9] = { 8.6370445090999749e+02, 0., 3.1603429399636133e+02, 
                      0., 8.6346178866552987e+02, 2.3044151161042186e+02, 
                      0., 0., 1.};
											*/
    
    /* chameleon
       double _cm[9] = { 833.61, 0, 351.11,
       0, 831.27, 228.38,
       0,  0,   1 };
    */
    
    Mat camMatrix = Mat(3,3,CV_64FC1,_cm);
    
    //          double _dc[] = {0.235, -0.878, 4.02, 0.013};
    //	double _dc[] = {0, 0, 0, 0,0};

    
    // firefly
    static double _dc[] = {-0.386530688588223, 0.160981787661065, 
		 											 0.00086827808373063, 0.00051791560238225, 
													  0};
    /*static double _dc[] = {-4.4608875275880239e-01, 1.4383997913318971e+00,
                           -1.1753032591690129e-04, 1.9647128820316859e-04,
                           -8.4978600160353448e+00 };
													 */
    /*	*/
    //      double _dc[] = { -1.5257215718379263e-01, -9.3288095506560065e-01,
    //                       7.0659389335978257e-04, 7.1318256767668705e-03,
    //                 4.8058340050954902e+00};
    
    
    //        cout << "O=" << O << endl;
    //        cout << "I=" << I << endl;
    //        cout << camMatrix << endl;
    
    cv::solvePnP(O, I, camMatrix, Mat(1,5,CV_64FC1,_dc), rvec, tvec, false);    
    for (int i = 0; i < 3; ++i) {
      q[i] = rvec.at<double>(i, 0);
      q[i+3] = tvec.at<double>(i, 0)/100;      
    }
    return true;
  }
  return false;
}
