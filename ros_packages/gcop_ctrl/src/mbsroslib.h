#ifndef MSROSLIB_H
#define MSROSLIB_H

#include "ros/ros.h"
#include <iomanip>
#include <iostream>
#include "gcop_comm/CtrlTraj.h"//msg for publishing ctrl trajectory
#include <gcop/urdf_parser.h>
#include "tf/transform_datatypes.h"
#include <gcop/se3.h>
#include "gcop/dmoc.h" //gcop dmoc header
#include "gcop/lqcost.h" //gcop lqr header
#include "gcop/mbscontroller.h"
#include "gcop_ctrl/MbsDMocInterfaceConfig.h"
#include <tf/transform_listener.h>
#include <XmlRpcValue.h>


using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<MbsState> MbsDmoc;//defining chaindmoc

class Mbstranslate 
{
	protected:
		//Pointer for mbs system
		boost::shared_ptr<Mbs> mbsmodel;
		//Pointer for Optimal Controller
		boost::shared_ptr<MbsDmoc> mbsdmoc;
		//MbsState final state
		boost::shared_ptr<MbsState> xf;
		//Cost lqcost
		boost::shared_ptr<LqCost<MbsState>> cost;
		boost::shared_ptr<MbsController> ctrl;
		int Nit = 1;//number of iterations for dmoc
		int N = 100;      // discrete trajectory segments
		string mbstype; // Type of system
		Matrix4d gposeroot_i; //inital inertial frame wrt the joint frame
	private:
		//functions
		void q2transform(geometry_msgs::Transform &transformmsg, Vector6d &bpose);
		void xml2vec(VectorXd &vec, XmlRpc::XmlRpcValue &my_list);
		void iterateCallback(const ros::TimerEvent & event);
		void paramreqcallback(gcop_ctrl::MbsDMocInterfaceConfig &config, uint32_t level);
	public:
		gcop_comm::CtrlTraj trajectory;
		MbsTranslate(string xmlstring);
};
#endif
