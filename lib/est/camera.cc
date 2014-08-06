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

/*
// firefly MTC
double _cm[9] = {623.77133283076, 0, 325.476798682001,
0, 624.840987895719, 180.290150726724,
0, 0, 1};

// firefly
double _dc[5] = {-0.386530688588223, 0.160981787661065, 
0.00086827808373063, 0.00051791560238225, 
0};
 */
	//static RNG rng(12345);

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


Camera::Camera(double bl, double imrows, double imcols,const double *K, const double *D):tolerance(2.0), kernel_size(2)
{
	this->bl = bl;
	//this->init_obj = false;
	filt = cv::Mat(imrows, imcols, CV_8UC1);
	rvec = cv::Mat::zeros(3,1, CV_64FC1);
	tvec = cv::Mat::zeros(3,1, CV_64FC1);
	tvec.at<float>(2,0) = 1;
	assert(K != 0);
	assert(D != 0);//Make sure there is some data
	double *camera_matrix = new double[9];
	double *distortion_vec = new double[5];
	memcpy(camera_matrix,K,9*sizeof(double));
	memcpy(distortion_vec,D,5*sizeof(double));
	camMatrix = Mat(3,3,CV_64FC1,camera_matrix);
	dccoeffs = Mat(1,5,CV_64FC1,distortion_vec);
	cout<<"CamMatrix: "<<endl<<camMatrix<<endl;
	cout<<"Distortion_Coeffs: "<<endl<<dccoeffs<<endl;
	pointfile.open("/home/gowtham/fuerte_workspace/sandbox/gcop/ros_packages/gcop_est/data/file1.dat",std::ofstream::out);
	cout<<"K "<<K[0]<<" "<<K[1]<<" "<<K[2]<<" "<<K[3]<<" "<<K[4]<<" "<<K[5]<<" "<<K[6]<<" "<<K[7]<<" "<<K[8]<<" "<<endl;
	//cout<<"D "<<D[0]<<" "<<D[1]<<" "<<D[2]<<" "<<D[3]<<" "<<D[4]<<endl;
}

Camera::~Camera()
{
	pointfile.close();
}



bool Camera::Pose(Vector6d &q, const cv::Mat &raw)//Input Raw cv image and output 6dof pose
{
	rvec = cv::Mat::zeros(3,1, CV_64FC1);
	tvec = cv::Mat::zeros(3,1, CV_64FC1);
	tvec.at<float>(2,0) = 1;
	threshold(raw, filt, 150, 255, THRESH_BINARY);
	blur( filt, filt, Size(kernel_size,kernel_size) );
	//medianBlur(filt,filt,3);//Median blur to remove salt and pepper noise

	vector< vector<Point> > contours;
	vector<Vec4i> hierarchy;
	findContours(filt, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);
	/*
	//Scalar color( 50, 50, 50 );
	cont =  Mat::zeros(filt.size(), CV_8UC3 );
	//raw.copyTo(cont);
	Scalar color( 50, 50, 50 );
	drawContours(cont, contours, -1, color,CV_FILLED);
	*/

	/*cout<<filt.size()<<endl;
	for( int i = 0; i< contours.size(); i++ )
	{
		Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
		drawContours( cont, contours, i, color, 2, 8, hierarchy, 0, Point());
	}
	*/

	/*
		 b a c

		 d   e
	 */

#define FLIP_AXIS

	//      float r = 4.5;
	//  float bl = 3.5;
#ifdef FLIP_AXIS


	static float _oprs[np][3] = {{(float)bl, 0, 0},
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


	Mat O(np, 3, CV_32FC1,_oprs);
	Mat I(np, 2, CV_32FC1);

	vector< vector<Point> >::iterator ci;
	if (contours.size() >= np) {
		if (contours.size() > np) {
			//cerr<<"More than 5 contours found"<<endl;//Warning if we have more than 5 contours
			//return false;
			//if(!init_obj)//If the object is not initialized do not get into the function unless we have exactly 5 contours
				//return false;
		}
		int i = 0;      
		int contnum = contours.size();
		vector<Point2f> center(contnum);
		vector<float> radius(contnum);
		//int pointindices[np];
		for(i=0;i<contnum;i++)
		{
			minEnclosingCircle(Mat(contours[i]),center[i],radius[i]);//Find the min enclosing circle
			vector<Point>::iterator pi;
			/*double x = 0, y = 0;
			for (pi = contours[i].begin(); pi != contours[i].end(); ++pi) {
				x += pi->x;
				y += pi->y;
			}
			x /= (contours[i]).size();
			y /= (contours[i]).size();
			*/
			//circle(cont, center[i], radius[i], Scalar( 255, 255, 0), 1);//Drawing circle
			//cout<<"Circle(x,y,r):\t"<<center[i].x<<"\t"<<center[i].y<<"\t"<<radius[i]<<"\t"<<x<<"\t"<<y<<endl;//The last two are for comparing
			//cout<<"Circle(x,y,r):\t"<<center[i].x<<"\t"<<center[i].y<<"\t"<<radius[i]<<endl;
			//center[i].x = x;
			//center[i].y = y;//Just to check how this will effect performance This did not improve the performance
		}
		//Now find the mode or most repeating radius (in our case it should be 5):
		//float tolerance = 1.5;
		int count3 = 0;
		//Algorithm for finding the closest 5 radii within tolerance (Not very optimal find a better algo TODO)
		for(int count1 = 0;count1 < contnum - 4; count1++)
		{
			//if(radius[count1] > 10 || radius[count1] < 2)
			//continue;//Just filter the radius to not search for radii which are too huge or too small
			iprs[count3].x = center[count1].x;
			iprs[count3].y = center[count1].y;//Add the first contour to the collection set before searching for parter contours
			count3++;
			//Go through all the radii count1+1 till end:
			for(int count2 = count1+1;count2 < contnum;count2++)
			{
				if(abs(radius[count2] - radius[count1]) < tolerance)
				{
					//cout<<"Counts for debug [count1, count2, count3, radius] :"<<count1<<"\t"<<count2<<"\t"<<count3<<"\t"<<radius[count1]<<endl;

					iprs[count3].x = center[count2].x;
					iprs[count3].y = center[count2].y;
					//pointindices[count3] = count2;
					count3+=1;//Incrementing count3
					if(count3 == np)
					{
						//Not handling this right now			cout<<"More than 5 contours with same radius"<<radius[pointindices[0]]<<endl;//Just printing one to see what radius it is
						break;//Done with 5 contours so no more searching
					}
				}
			}
			if(count3 == np)
				break;//Done searching
			else
				count3 = 0;//Else start recording points again
		}
		if(count3 < np)//If we found the 5 contours then use them
		{
			return false;//Not enough contours
		}
		/*
		else
		{
			cout<<"Found 5 contours Hurray !"<<endl;
		}
		*/
		//This does not handle the case of more than 5 points with similar radius and hopefull we will not get many of such cases
		//if(!init_obj)
		//init_obj = true;
		//cout<<"iprs: \n"<<iprs[0].x<<"\t"<<iprs[0].y<<"\n"<<iprs[1].x<<"\t"<<iprs[1].y<<"\n"<<iprs[2].x<<"\t"<<iprs[2].y<<"\n"<<iprs[3].x<<"\t"<<iprs[3].y<<"\n"<<iprs[4].x<<"\t"<<iprs[4].y<<"\n"<<endl;
		double min = numeric_limits<double>::max();
		//Reorder the points
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


		//          double _dc[] = {0.235, -0.878, 4.02, 0.013};
		//	double _dc[] = {0, 0, 0, 0,0};


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

		cv::solvePnP(O, I, camMatrix, dccoeffs, rvec, tvec, false);    
		//Write the 5 points to pointfile:
		for(int i = 0;i<np;i++)
		{
			pointfile<<iprs[i].x<<"\t"<<iprs[i].y<<"\t";
		}
		pointfile<<endl;
		for (int i = 0; i < 3; ++i) {
			q[i] = rvec.at<double>(i, 0);
			q[i+3] = tvec.at<double>(i, 0)/100;      
		}
		return true;
	}
	//init_obj = false;//The object is reinitialized everytime we loose it
	return false;
}
bool Camera::Pose(Vector6d &q, const cv::Mat &raw, cv::Mat &cont)//Input Raw cv image and output 6dof pose
{
	rvec = cv::Mat::zeros(3,1, CV_64FC1);
	tvec = cv::Mat::zeros(3,1, CV_64FC1);
	tvec.at<float>(2,0) = 1;
	threshold(raw, filt, 150, 255, THRESH_BINARY);
	blur( filt, filt, Size(kernel_size,kernel_size) );
	//medianBlur(filt,filt,3);//Median blur to remove salt and pepper noise

	vector< vector<Point> > contours;
	vector<Vec4i> hierarchy;
	findContours(filt, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);
	//Scalar color( 50, 50, 50 );
	cont =  Mat::zeros(filt.size(), CV_8UC3 );
	//raw.copyTo(cont);
	Scalar color( 50, 50, 50 );
	drawContours(cont, contours, -1, color,CV_FILLED);

	/*cout<<filt.size()<<endl;
	for( int i = 0; i< contours.size(); i++ )
	{
		Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
		drawContours( cont, contours, i, color, 2, 8, hierarchy, 0, Point());
	}
	*/

	/*
		 b a c

		 d   e
	 */

#define FLIP_AXIS

	//      float r = 4.5;
	//  float bl = 3.5;
#ifdef FLIP_AXIS


	static float _oprs[np][3] = {{(float)bl, 0, 0},
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


	Mat O(np, 3, CV_32FC1,_oprs);
	Mat I(np, 2, CV_32FC1);

	vector< vector<Point> >::iterator ci;
	if (contours.size() >= np) {
		if (contours.size() > np) {
			//cerr<<"More than 5 contours found"<<endl;//Warning if we have more than 5 contours
			//return false;
			//if(!init_obj)//If the object is not initialized do not get into the function unless we have exactly 5 contours
				//return false;
		}
		int i = 0;      
		int contnum = contours.size();
		vector<Point2f> center(contnum);
		vector<float> radius(contnum);
		//int pointindices[np];
		for(i=0;i<contnum;i++)
		{
			minEnclosingCircle(Mat(contours[i]),center[i],radius[i]);//Find the min enclosing circle
			vector<Point>::iterator pi;
			double x = 0, y = 0;
			/*for (pi = contours[i].begin(); pi != contours[i].end(); ++pi) {
				x += pi->x;
				y += pi->y;
			}
			x /= (contours[i]).size();
			y /= (contours[i]).size();
			*/
			circle(cont, center[i], radius[i], Scalar( 255, 255, 0), 1);//Drawing circle
			//cout<<"Circle(x,y,r):\t"<<center[i].x<<"\t"<<center[i].y<<"\t"<<radius[i]<<"\t"<<x<<"\t"<<y<<endl;//The last two are for comparing
			//cout<<"Circle(x,y,r):\t"<<center[i].x<<"\t"<<center[i].y<<"\t"<<radius[i]<<endl;
			//center[i].x = x;
			//center[i].y = y;//Just to check how this will effect performance This did not improve the performance
		}
		//Now find the mode or most repeating radius (in our case it should be 5):
		//float tolerance = 1.5;
		int count3 = 0;
		//Algorithm for finding the closest 5 radii within tolerance (Not very optimal find a better algo TODO)
		for(int count1 = 0;count1 < contnum - 4; count1++)
		{
			//if(radius[count1] > 10 || radius[count1] < 2)
			//continue;//Just filter the radius to not search for radii which are too huge or too small
			iprs[count3].x = center[count1].x;
			iprs[count3].y = center[count1].y;//Add the first contour to the collection set before searching for parter contours
			count3++;
			//Go through all the radii count1+1 till end:
			for(int count2 = count1+1;count2 < contnum;count2++)
			{
				if(abs(radius[count2] - radius[count1]) < tolerance)
				{
					//cout<<"Counts for debug [count1, count2, count3, radius] :"<<count1<<"\t"<<count2<<"\t"<<count3<<"\t"<<radius[count1]<<endl;

					iprs[count3].x = center[count2].x;
					iprs[count3].y = center[count2].y;
					//pointindices[count3] = count2;
					count3+=1;//Incrementing count3
					if(count3 == np)
					{
						//Not handling this right now			cout<<"More than 5 contours with same radius"<<radius[pointindices[0]]<<endl;//Just printing one to see what radius it is
						break;//Done with 5 contours so no more searching
					}
				}
			}
			if(count3 == np)
				break;//Done searching
			else
				count3 = 0;//Else start recording points again
		}
		if(count3 < np)//If we found the 5 contours then use them
		{
			return false;//Not enough contours
		}
		/*
		else
		{
			cout<<"Found 5 contours Hurray !"<<endl;
		}
		*/
		//This does not handle the case of more than 5 points with similar radius and hopefull we will not get many of such cases
		/*
			 for (ci = contours.begin(); ci != contours.end(); ++ci) { 
			 vector<Point>::iterator pi;
			 double x = 0, y = 0, xsqr = 0, ysqr = 0;
			 for (pi = ci->begin(); pi != ci->end(); ++pi) {
			 x += pi->x;
			 y += pi->y;
		//xsqr += (pi->x)*(pi->x);
		//ysqr += (pi->y)*(pi->y);
		}
		x /= (*ci).size();
		y /= (*ci).size();
		//Fitting minenclosing circle to the contours:
		minEnclosingCircle(Mat(contours[i]),center,radius);//Find the min enclosing circle
		//RotatedRect bounding_rect = minAreaRect( Mat(ci));

		//xsqr /= (*ci).size();
		//ysqr /= (*ci).size();
		//double radius = sqrt((1.0f/(*ci).size())*(xsqr - x*x + ysqr - y*y));
		cout<<"X Y:  "<<x<<"\t"<<y<<endl;
		cout<<"Circle(x,y,r):\t"<<center.x<<"\t"<<center.y<<"\t"<<radius<<endl;
		if(i < np )//Have to deal with small dot contours greater than 5
		{
		iprs[i].x = x;
		iprs[i].y = y;
		}
		++i;
		 *** previous comment begin	
		 int visitedone = 0;
		 double mindist = numeric_limits<double>::max();
		 for(int count1 = 0;count1 < np;count1++)
		 {
		 if(!visited[count1])
		 {
		 double distance = sqrt((iprs[count1].x - x)*(iprs[count1].x - x) + (iprs[count1].y - y)*(iprs[count1].y - y));
		 visitedone = distance < mindist?count1:visitedone;
		 mindist = distance < mindist?distance:mindist;
		 }
		 }
		 cout<<mindist<<endl;
		 if(mindist < 100 || !init_obj)
		 {
		 visited[visitedone] = true;//that point is visited 
		 cout<<" Visited one: "<<visitedone<<endl;
		 iprs[visitedone].x = x;
		 iprs[visitedone].y = y;
		 }
		 Check areas of the contours thats the best as of now
		 *** previous commented end
		//check if all are visited TODO
		//if (i < np) {
		//		cout<<iprs[i].x<<"\t"<<iprs[i].y<<endl;
		//cout<<i<<" iprs: \n"<<iprs[0].x<<"\t"<<iprs[0].y<<"\n"<<iprs[1].x<<"\t"<<iprs[1].y<<"\n"<<iprs[2].x<<"\t"<<iprs[2].y<<"\n"<<iprs[3].x<<"\t"<<iprs[3].y<<"\n"<<iprs[4].x<<"\t"<<iprs[4].y<<"\n"<<endl;
		//cout<<"Distance["<<i<<"]: "<<sqrt((iprs[i].x - x)*(iprs[i].x - x) + (iprs[i].y - y)*(iprs[i].y - y))<<endl;
		//	iprs[i].x = x;
		//	iprs[i].y = y;
		//	++i;
		//}
		}
		if(i != 5)
		{
		cout<<i<<endl;
		return false;//This means it cannot find 5 contours which are of the right radius
		}
		 */
		//if(!init_obj)
		//init_obj = true;
		//cout<<"iprs: \n"<<iprs[0].x<<"\t"<<iprs[0].y<<"\n"<<iprs[1].x<<"\t"<<iprs[1].y<<"\n"<<iprs[2].x<<"\t"<<iprs[2].y<<"\n"<<iprs[3].x<<"\t"<<iprs[3].y<<"\n"<<iprs[4].x<<"\t"<<iprs[4].y<<"\n"<<endl;
		double min = numeric_limits<double>::max();
		//Reorder the points
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


		//          double _dc[] = {0.235, -0.878, 4.02, 0.013};
		//	double _dc[] = {0, 0, 0, 0,0};


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

		cv::solvePnP(O, I, camMatrix, dccoeffs, rvec, tvec, false);    
		//Write the 5 points to pointfile:
		for(int i = 0;i<np;i++)
		{
			pointfile<<iprs[i].x<<"\t"<<iprs[i].y<<"\t";
		}
		pointfile<<endl;
		for (int i = 0; i < 3; ++i) {
			q[i] = rvec.at<double>(i, 0);
			q[i+3] = tvec.at<double>(i, 0)/100;      
		}
		return true;
	}
	//init_obj = false;//The object is reinitialized everytime we loose it
	return false;
}
