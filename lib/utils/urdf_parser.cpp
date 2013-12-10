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

#include <boost/algorithm/string.hpp>
#include <vector>
#include "urdf_parser.h"
#include "se3.h"

using namespace std;
using namespace Eigen;

namespace gcop_urdf{
/*
gcop::Matrix4d diffpose(Pose &posej_p,Pose &posei_p)
{
	//p parent i inertial j joint we need j_i joint wrt i = j_i = j_p * p_i = j_p *(i_p)^-1
	gcop::Matrix4d gposei_p;
	gcop::Matrix4d gposep_i;
	gcop::Matrix4d gposej_p;
	gcop::Matrix4d gposej_i;
	double r, p, y;
	posej_p.rotation.getRPY(r,p,y);
	gcop::SE3::Instance().rpyxyz2g(gposej_p, Vector3d(r,p,y), Vector3d(posej_p.position.x,posej_p.position.y,posej_p.position.z));
	//get i_p
	posei_p.rotation.getRPY(r,p,y);
	gcop::SE3::Instance().rpyxyz2g(gposei_p, Vector3d(r,p,y), Vector3d(posei_p.position.x,posei_p.position.y,posei_p.position.z));
	gcop::SE3::Instance().inv(gposep_i,gposei_p);

	gposej_i = gposej_p*gposep_i;

	return gposej_i;
}

//link is the parent link and level is the index of the link in pis;
//bodies, pis and joints are the list to be used in the mbs class
void walkTree(boost::shared_ptr<const Link> link, int level,int &index,boost::shared_ptr<gcop::Mbs> mbs)
{
	//int count = 0;
	for (vector<boost::shared_ptr<Link> >::const_iterator child = link->child_links.begin(); child != link->child_links.end(); child++)
	{
		if (*child)
		{
			boost::shared_ptr<const Link> childlink = *child;
			// first grandchild

			//get parent joint of the child:
			//boost::shared_ptr<Joint> parentj = child->
			//create joint
			//gcop::Joint gcopjoint;
			//based on type  of joint decide the axis
		 switch(childlink->parent_joint.type)
		 {
			 case Joint::REVOLUTE :
				 mbs->joints[index].a.setZero();
				 mbs->joints[index].a[0]  = childlink->parent_joint->axis.x;
				 mbs->joints[index].a[1]  = childlink->parent_joint->axis.y;
				 mbs->joints[index].a[2]  = childlink->parent_joint->axis.z;
				 break;
			 case Joint::PRISMATIC :
				 mbs->joints[index].a.setZero();
				 mbs->joints[index].a[3]  = childlink->parent_joint->axis.x;
				 mbs->joints[index].a[4]  = childlink->parent_joint->axis.y;
				 mbs->joints[index].a[5]  = childlink->parent_joint->axis.z;
				 break;
		 }
		 //decide the transformations:
			mbs->joints[index].gp = diffpose(childlink->parent_joint->parent_to_joint_origin_transform,link->inertial->origin);
			double r,p,y;
			childlink->inertial->origin.rotation.getRPY(r,p,y);

			gcop::SE3::Instance().rpyxyz2g(mbs->joints[index].gc,Vector3d(r,p,y),Vector3d(childlink->inertial->origin.position.x,childlink->inertial->origin.position.y,childlink->inertial->origin.position.z));

			cout<<"level: "<<level<<endl;
			
			index++;

			//create body3d from the current link and make a joint connecting them. This joint is the parent joint of the child link
			
			//gcop::Body3d<> body;
			//copy m
			mbs->links[index].m = childlink->inertial->mass;
			//copy J
			mbs->links[index].J(0) = childlink->inertial->ixx;
			mbs->links[index].J(1) = childlink->inertial->iyy;
			mbs->links[index].J(2) = childlink->inertial->izz;//we assume the axes are aligned with principal axes
			// make ds = 0;
			mbs->links[index].ds = Vector3d(0,0,0);

			mbs->pis[index] = level; //added parent index to the list.


			walkTree(*child,index,index,mbs);
		}
		else
		{
			for(int j=0;j<level;j++) cout << " "; //indent
			cout << "root link: " << link->name << " has a null child!" << *child << endl;
		}
	}
}
*/
boost::shared_ptr<gcop::Mbs> mbsgenerator(const string &xml_string)
{
	//parse the xml file into a urdf model
	boost::shared_ptr<ModelInterface> urdfmodel = parseURDF(xml_string);

	if (!urdfmodel){
		cerr << "ERROR: Model Parsing the xml failed" << endl;
		return boost::shared_ptr<gcop::Mbs>();
	}
	cout << "robot name is: " << urdfmodel->getName() << endl;

	// get root link
	boost::shared_ptr<const Link> root_link = urdfmodel->getRoot();
	if (!root_link) return boost::shared_ptr<gcop::Mbs>();

	boost::shared_ptr<Link> baselink = root_link->child_links[0];
	if(!root_link) return boost::shared_ptr<gcop::Mbs>();

	cout<<urdfmodel->links_.size() <<"\t"<<urdfmodel->joints_.size()<<endl;


	//create Mbs with the number of links and joints
	boost::shared_ptr<gcop::Mbs> mbs (new gcop::Mbs(urdfmodel->links_.size(),6 + (urdfmodel->joints_.size()-1)));//hack for now

	int index = 0;

	//add the first link to mbs
	mbs->links[index].m = baselink->inertial->mass;
	//copy J
	mbs->links[index].J(0) = baselink->inertial->ixx;
	mbs->links[index].J(1) = baselink->inertial->iyy;
	mbs->links[index].J(2) = baselink->inertial->izz;//we assume the axes are aligned with principal axes
	// make ds = 0;
	mbs->links[index].ds = Vector3d(0,0,0);

	mbs->pis[index] = -1; //added parent index to the list.


	//walkTree(root_link,-1,index,mbs);
	cout<<index<<endl;
	cout<<"Number of links parsed: "<<bodies.size()<<endl;
	cout<<" Number of joints parsed: "<<joints.size()<<endl;
	//cout<<"Pis: "<<pis<<endl;
	return mbs;
}



boost::shared_ptr<ModelInterface>  parseURDF(const string &xml_string)
{
  boost::shared_ptr<ModelInterface> model(new ModelInterface);
  model->clear();

  TiXmlDocument xml_doc;
  xml_doc.Parse(xml_string.c_str());

  TiXmlElement *robot_xml = xml_doc.FirstChildElement("robot");
  if (!robot_xml)
  {
    printf("Could not find the 'robot' element in the xml file");
    model.reset();
    return model;
  }

  // Get robot name
  const char *name = robot_xml->Attribute("name");
  if (!name)
  {
    printf("No name given for the robot.");
    model.reset();
    return model;
  }
  model->name_ = string(name);

  // Get all Material elements
  for (TiXmlElement* material_xml = robot_xml->FirstChildElement("material"); material_xml; material_xml = material_xml->NextSiblingElement("material"))
  {
    boost::shared_ptr<Material> material;
    material.reset(new Material);

    if (material->initXml(material_xml))
    {
      if (model->getMaterial(material->name))
      {
        printf("material '%s' is not unique.", material->name.c_str());
        material.reset();
        model.reset();
        return model;
      }
      else
      {
        model->materials_.insert(make_pair(material->name,material));
        printf("successfully added a new material '%s'", material->name.c_str());
      }
    }
    else
    {
      printf("material xml is not initialized correctly");
      material.reset();
      model.reset();
      return model;
    }
  }

  // Get all Link elements
  for (TiXmlElement* link_xml = robot_xml->FirstChildElement("link"); link_xml; link_xml = link_xml->NextSiblingElement("link"))
  {
    boost::shared_ptr<Link> link;
    link.reset(new Link);

    if (link->initXml(link_xml))
    {
      if (model->getLink(link->name))
      {
        printf("link '%s' is not unique.", link->name.c_str());
        model.reset();
        return model;
      }
      else
      {
        // set link visual material
        printf("setting link '%s' material", link->name.c_str());
        if (link->visual)
        {
          if (!link->visual->material_name.empty())
          {
            if (model->getMaterial(link->visual->material_name))
            {
              printf("setting link '%s' material to '%s'", link->name.c_str(),link->visual->material_name.c_str());
              link->visual->material = model->getMaterial( link->visual->material_name.c_str() );
            }
            else
            {
              if (link->visual->material)
              {
                printf("link '%s' material '%s' defined in Visual.", link->name.c_str(),link->visual->material_name.c_str());
                model->materials_.insert(make_pair(link->visual->material->name,link->visual->material));
              }
              else
              {
                printf("link '%s' material '%s' undefined.", link->name.c_str(),link->visual->material_name.c_str());
                model.reset();
                return model;
              }
            }
          }
        }

        model->links_.insert(make_pair(link->name,link));
        printf("successfully added a new link '%s'", link->name.c_str());
      }
    }
    else
    {
      printf("link xml is not initialized correctly");
      model.reset();
      return model;
    }
  }
  if (model->links_.empty()){
    printf("No link elements found in urdf file");
    model.reset();
    return model;
  }

  // Get all Joint elements
  for (TiXmlElement* joint_xml = robot_xml->FirstChildElement("joint"); joint_xml; joint_xml = joint_xml->NextSiblingElement("joint"))
  {
    boost::shared_ptr<Joint> joint;
    joint.reset(new Joint);

    if (joint->initXml(joint_xml))
    {
      if (model->getJoint(joint->name))
      {
        printf("joint '%s' is not unique.", joint->name.c_str());
        model.reset();
        return model;
      }
      else
      {
        model->joints_.insert(make_pair(joint->name,joint));
        printf("successfully added a new joint '%s'", joint->name.c_str());
      }
    }
    else
    {
      printf("joint xml is not initialized correctly");
      model.reset();
      return model;
    }
  }


  // every link has children links and joints, but no parents, so we create a
  // local convenience data structure for keeping child->parent relations
  map<string, string> parent_link_tree;
  parent_link_tree.clear();

  // building tree: name mapping
  if (!model->initTree(parent_link_tree))
  {
    printf("failed to build tree");
    model.reset();
    return model;
  }

  // find the root link
  if (!model->initRoot(parent_link_tree))
  {
    printf("failed to find root link");
    model.reset();
    return model;
  }
 
  return model;
}

}

