/*********************************************************************
* Software License Agreement (BSD License)
* 
*  Copyright (c) 2008, Willow Garage, Inc.
*  All rights reserved.
* 
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
* 
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Willow Garage nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
* 
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Wim Meeussen */
/* ReAuthoring: Gowtham Garimella */

#include <boost/algorithm/string.hpp>
#include <ros/console.h>
#include <vector>
#include "urdf_parser/urdf_parser.h"


namespace urdf{



boost::shared_ptr<ModelInterface>  parseURDF(const std::string &xml_string)
{
  //boost::shared_ptr<ModelInterface> model(new ModelInterface);
  boost::shared_ptr<gcop::Mbs> model(new gcop::Mbs());//do this after we know number of joints and links
	std::vector<gcop::Body3d> linklist;
  //model->clear();
	/*Create Multi Body system and return that*/
	/**/
	/**/

  TiXmlDocument xml_doc;
  xml_doc.Parse(xml_string.c_str());

  TiXmlElement *robot_xml = xml_doc.FirstChildElement("robot");
  if (!robot_xml)
  {
    ROS_ERROR("Could not find the 'robot' element in the xml file");
    model.reset();
    return model;
  }

 	/*Convert these link elements to mbs Body3d */
  // Get all Link elements
  for (TiXmlElement* link_xml = robot_xml->FirstChildElement("link"); link_xml; link_xml = link_xml->NextSiblingElement("link"))
  {
    boost::shared_ptr<Link> link;
		gcop::Body3d gcoplink;

		
    link.reset(new Link);

    if (link->initXml(link_xml))
    {
      if (model->getLink(link->name))
      {
        ROS_ERROR("link '%s' is not unique.", link->name.c_str());
        model.reset();
        return model;
      }
      else
      {
        if (link->inertial)
				{

				}
        model->links_.insert(make_pair(link->name,link));
        ROS_DEBUG("successfully added a new link '%s'", link->name.c_str());
      }
    }
    else
    {
      ROS_ERROR("link xml is not initialized correctly");
      //model.reset();
      //return model;
    }
  }
  if (model->links_.empty()){
    ROS_ERROR("No link elements found in urdf file");
    model.reset();
    return model;
  }

	/* Convert all these joints to mbs joints */
  // Get all Joint elements
  for (TiXmlElement* joint_xml = robot_xml->FirstChildElement("joint"); joint_xml; joint_xml = joint_xml->NextSiblingElement("joint"))
  {
    boost::shared_ptr<Joint> joint;
    joint.reset(new Joint);

    if (joint->initXml(joint_xml))
    {
      if (model->getJoint(joint->name))
      {
        ROS_ERROR("joint '%s' is not unique.", joint->name.c_str());
        model.reset();
        return model;
      }
      else
      {
        model->joints_.insert(make_pair(joint->name,joint));
        ROS_DEBUG("successfully added a new joint '%s'", joint->name.c_str());
      }
    }
    else
    {
      ROS_ERROR("joint xml is not initialized correctly");
      model.reset();
      return model;
    }
  }


  // every link has children links and joints, but no parents, so we create a
  // local convenience data structure for keeping child->parent relations
  std::map<std::string, std::string> parent_link_tree;
  parent_link_tree.clear();

  // building tree: name mapping
  if (!model->initTree(parent_link_tree))
  {
    ROS_ERROR("failed to build tree");
    model.reset();
    return model;
  }

  // find the root link
  if (!model->initRoot(parent_link_tree))
  {
    ROS_ERROR("failed to find root link");
    model.reset();
    return model;
  }
 
  return model;
}

}

