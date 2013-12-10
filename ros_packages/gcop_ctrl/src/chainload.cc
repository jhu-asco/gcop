#include "ros/ros.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <dynamic_reconfigure/server.h>
#include <gcop/urdf_parser.h>

using namespace std;
using namespace Eigen;
using namespace gcop;

int main(int argc, char** argv)
{
	ros::init(argc, argv, "chainload");
	ros::NodeHandle n;
	//get parameter for xml_string:
	string xml_string, xml_filename;
	if(!ros::param::get("/robot_description", xml_filename))
	{
		ROS_ERROR("Could not fetch xml file parameter");
		return 0;
	}
	cout<<xml_filename<<endl;
	fstream xml_file(xml_filename.c_str(), fstream::in);
	while ( xml_file.good() )
	{
		string line;
		getline( xml_file, line);
		xml_string += (line + "\n");
	}
	xml_file.close();

	boost::shared_ptr<Mbs> mbsmodel = gcop_urdf::mbsgenerator(xml_string);

	//check if mbsmodel is good by using it.



	//boost::shared_ptr<gcop_urdf::ModelInterface> robot = gcop_urdf::parseURDF(xml_string);
	/*
	if (!robot){
		cerr << "ERROR: Model Parsing the xml failed" << endl;
		return -1;
	}
	cout << "robot name is: " << robot->getName() << endl;

	// get info from parser
	cout << "---------- Successfully Parsed XML ---------------" << endl;
	// get root link
	boost::shared_ptr<const Link> root_link=robot->getRoot();
	if (!root_link) return -1;

	cout << "root Link: " << root_link->name << " has " << root_link->child_links.size() << " child(ren)" << endl;
	*/




	return 0;
}
