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
#include <Eigen/Eigenvalues>
#include "se3.h"

using namespace std;
using namespace Eigen;

namespace gcop_urdf{
gcop::Matrix4d diffpose(Pose &posej_p,Pose &posei_p)
{
	//p parent i inertial j joint we need j_i joint wrt i = j_i = j_p * p_i = j_p *(i_p)^-1
	/*gcop::Matrix4d gposei_p;
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
	*/
	Pose posej_i = (posei_p.GetInverse())*posej_p;
	gcop::Matrix4d gposej_i;
	gcop::Vector7d posevec;
	posevec<<posej_i.rotation.w,posej_i.rotation.x,posej_i.rotation.y,posej_i.rotation.z,posej_i.position.x,posej_i.position.y,posej_i.position.z;
	gcop::SE3::Instance().quatxyz2g(gposej_i,posevec);
	return gposej_i;
}

void transformtoprincipal(boost::shared_ptr<Link> link)
{
	Vector3d offdiag(link->inertial->ixy,link->inertial->ixz,link->inertial->iyz);
	if(offdiag.norm() < 1e-10)
		return;
	gcop::Matrix3d linkinertia;
	linkinertia(0,0) = link->inertial->ixx; linkinertia(0,1) = link->inertial->ixy; linkinertia(0,2) = link->inertial->ixz; linkinertia(1,1) = link->inertial->iyy;linkinertia(1,2) = link->inertial->iyz;linkinertia(2,2) = link->inertial->izz;
	//solve for principal inertial vector:
	SelfAdjointEigenSolver<Matrix3d> es(linkinertia.selfadjointView<Eigen::Upper>());
	Vector3d finalprincipal = es.eigenvalues();
	cout<<"Final principal inertia values"<<endl<<es.eigenvalues()<<endl;
	link->inertial->ixx = finalprincipal(0);	link->inertial->ixy = 0; link->inertial->ixz = 0; link->inertial->iyy = finalprincipal(1); link->inertial->iyz = 0; link->inertial->izz = finalprincipal(2);	
	//add the new transform  for new COM 
	Matrix3d principaltransform = es.eigenvectors();
	cout<<"Principal transform"<<endl<<principaltransform<<endl;
	if((principaltransform.determinant()) < 0)
	{
		principaltransform.col(0).swap(principaltransform.col(1));
		link->inertial->ixx = finalprincipal(1);	link->inertial->ixy = 0; link->inertial->ixz = 0; link->inertial->iyy = finalprincipal(0); link->inertial->iyz = 0; link->inertial->izz = finalprincipal(2);	
		cout<<"Final principal inertia values "<<finalprincipal(1)<<" "<<finalprincipal(0)<<" "<<finalprincipal(2)<<endl;
		cout<<"principaltransform "<<endl<<principaltransform<<endl;
	}
	cout<<"Determinant of evectors: "<<principaltransform.determinant()<<endl;
	//Matrix3d R = principaltransform.topLeftCorner<3,3>();
	Matrix3d zero = principaltransform*principaltransform.transpose() - Matrix3d::Identity();
	//cout<<"trial lpnorm "<<zero.determinant()<<endl;
	//assert statement for checking rotation matrix:
	assert(abs(zero.determinant()) < 1e-3);

	gcop::Vector4d vrotfi_pi;
	gcop::SO3::Instance().g2quat(vrotfi_pi,principaltransform);
	Rotation rotfi_pi;
	//Pose posefi_pi;
	rotfi_pi.w = vrotfi_pi(0); rotfi_pi.x = vrotfi_pi(1); rotfi_pi.y = vrotfi_pi(2); rotfi_pi.z = vrotfi_pi(3); 
	//cout<<"posefi_pi "<<posefi_pi<<endl;
	//do inverse conversion for checking integrity of g2quat:
	//	Matrix4d ptransform;
	//	gcop::SE3::Instance().quatxyz2g(ptransform,vposefi_pi);
	//cout<<"ptransform "<<ptransform<<endl;

	link->inertial->origin.rotation = link->inertial->origin.rotation*rotfi_pi;
}

//Convenience function for one unit of combining parent link and child link
void combineinertia(boost::shared_ptr<Link> clink, boost::shared_ptr<Link> plink,Pose posec_p)
{
	Pose poseci_p =  posec_p*(clink->inertial->origin);
	Pose poseci_pi = plink->inertial->origin.GetInverse()*poseci_p;
	cout<<"Poseci_p"<<poseci_p<<endl;
	cout<<"Poseci_pi"<<poseci_pi<<endl;


	//	Inertial finalinertia;//implicit copy constructor is enough here
	//finalinertia.position = (pmass*plink->inertial->origin.position + cmass*poseci_p.position)/(pmass + cmass);

	gcop::Matrix4d gposeci_pi; 
	gcop::Vector7d posevec;
	posevec<<poseci_pi.rotation.w,poseci_pi.rotation.x,poseci_pi.rotation.y,poseci_pi.rotation.z,poseci_pi.position.x,poseci_pi.position.y,poseci_pi.position.z;
	gcop::SE3::Instance().quatxyz2g(gposeci_pi,posevec);
	//convert inertia of childlink to pi frame
	gcop::Matrix3d rotationci_pi = gposeci_pi.topLeftCorner<3,3>();
	gcop::Matrix3d childinertia;
	cout<<"rotationci_pi"<<endl<<rotationci_pi<<endl;
	childinertia(0,0) = clink->inertial->ixx; childinertia(0,1) = clink->inertial->ixy; childinertia(0,2) = clink->inertial->ixz; childinertia(1,1) = clink->inertial->iyy;childinertia(1,2) = clink->inertial->iyz;childinertia(2,2) = clink->inertial->izz;

	cout<<"Childinertia before rotation"<<endl<<childinertia<<endl;

	childinertia.triangularView<Upper>()  = rotationci_pi * childinertia.selfadjointView<Eigen::Upper>() *rotationci_pi.transpose();

	cout<<"Childinertia after rotation"<<endl<<childinertia<<endl;


	gcop::Matrix3d parentinertia;
	parentinertia(0,0) = plink->inertial->ixx; parentinertia(0,1) = plink->inertial->ixy; parentinertia(0,2) = plink->inertial->ixz; parentinertia(1,1) = plink->inertial->iyy;parentinertia(1,2) = plink->inertial->iyz;parentinertia(2,2) = plink->inertial->izz;
	cout<<"ParentInertia "<<endl<<parentinertia<<endl;

	//Using parallel axis theorem:
	double pmass = plink->inertial->mass, cmass = clink->inertial->mass;
	cout<<"pmass "<<pmass<<" cmass "<<cmass<<endl;
	Vector3d vecfi_pi = gposeci_pi.topRightCorner<3,1>()*(pmass/(cmass + pmass)); 
	Vector3d vecfi_ci = vecfi_pi - gposeci_pi.topRightCorner<3,1>();

	cout<<"vecfi_pi "<<vecfi_pi.transpose()<<" vecfi_ci "<<vecfi_ci.transpose()<<endl;

	gcop::Matrix3d finalinertia;
	finalinertia.triangularView<Upper>() = parentinertia + pmass*(vecfi_pi.squaredNorm()*Matrix3d::Identity() - vecfi_pi*vecfi_pi.transpose())+ childinertia + cmass*(vecfi_ci.squaredNorm()*Matrix3d::Identity() - vecfi_ci*vecfi_ci.transpose()) ;

	//assign the final inertia values to the parentlink
	plink->inertial->mass = pmass + cmass;
	plink->inertial->ixx = finalinertia(0,0);	plink->inertial->ixy = finalinertia(0,1); plink->inertial->ixz = finalinertia(0,2); plink->inertial->iyy = finalinertia(1,1); plink->inertial->iyz = finalinertia(1,2); plink->inertial->izz = finalinertia(2,2);	
	Vector3 pvecfi_pi(vecfi_pi(0),vecfi_pi(1),vecfi_pi(2));
	plink->inertial->origin.position = plink->inertial->origin.position + plink->inertial->origin.rotation*pvecfi_pi;//put the frame at the center of mass
}
//Assign the inertial origin of the parent link to all the child links. This needs to be done after aggregation to ensure that all the fixed child links have access to the latest inertial information
void assign(boost::shared_ptr<const Link> link, boost::shared_ptr<Link> parentlink,Pose cumpose)
{
	for (vector<boost::shared_ptr<Link> >::const_iterator child = link->child_links.begin(); child != link->child_links.end(); child++)
	{
		if (*child)
		{
			boost::shared_ptr<Link> childlink = *child;
			// first grandchild

			//get parent joint of the child:
			//based on type  of joint decide to combine or just return 
			if(childlink->parent_joint->type != Joint::FIXED)
			{
				continue;// if joint is not fixed then dont assign 
			}
			//Now we have a fixed joint try assigning the inertial frame 
			Pose posec_p = cumpose*childlink->parent_joint->parent_to_joint_origin_transform;
			childlink->inertial->origin = posec_p.GetInverse()*parentlink->inertial->origin;
			cout<<"Posec_p "<<posec_p<<endl;

			gcop::Matrix4d gposei_c;
			gcop::Vector7d posevec;
			posevec<<childlink->inertial->origin.rotation.w,childlink->inertial->origin.rotation.x,childlink->inertial->origin.rotation.y,childlink->inertial->origin.rotation.z,childlink->inertial->origin.position.x,childlink->inertial->origin.position.y,childlink->inertial->origin.position.z;
			gcop::SE3::Instance().quatxyz2g(gposei_c,posevec);
			cout<<"gposei_c"<<endl<<gposei_c<<endl;
			cout<<"Posei_c "<<posevec.transpose()<<endl;
			childlink->visited = true;
			assign(childlink,parentlink,posec_p);
		}
		else
		{
			cout << "root link: " << link->name << " has a null child!" << *child << endl;
		}
	}
}
//recursive combining algorithm for fixed joints. parentlink is the final resulting link to wich all the links are merged
void aggregate(boost::shared_ptr<const Link> link, boost::shared_ptr<Link> parentlink,Pose cumpose)
{
for (vector<boost::shared_ptr<Link> >::const_iterator child = link->child_links.begin(); child != link->child_links.end(); child++)
	{
		if (*child)
		{
			boost::shared_ptr<Link> childlink = *child;
			// first grandchild

	 	 //get parent joint of the child:
		 //based on type  of joint decide to combine or just return 
		 if(childlink->parent_joint->type != Joint::FIXED)
		 {
			 continue;// if joint is not fixed then dont combine
		 }
		 cout<<childlink->name<<endl;
		 //Now we have a fixed joint try combining them 
		 Pose cumposenext = cumpose*childlink->parent_joint->parent_to_joint_origin_transform;
		 combineinertia(childlink,parentlink,cumposenext);
		 //then aggregate the fixed links of childlink:
		 aggregate(childlink,parentlink,cumposenext);
		}
		else
		{
				cout << "root link: " << link->name << " has a null child!" << *child << endl;
		}
	}
}

//link is the parent link and level is the index of the link in pis;
//bodies, pis and joints are the list to be used in the mbs class
void walkTree(boost::shared_ptr<Link> link, int level,int &index,boost::shared_ptr<gcop::Mbs> mbs)
{	
	//int count = 0;
	for (vector<boost::shared_ptr<Link> >::const_iterator child = link->child_links.begin(); child != link->child_links.end(); child++)
	{
		if (*child)
		{
			boost::shared_ptr<Link> childlink = *child;
			// first grandchild
			if(!childlink->visited)
			{
				childlink->visited = true;//marks the link as visited for aggregation purposes
				Pose idpose;
				idpose.clear();
				aggregate(childlink,childlink,idpose);//aggregates all the fixed joint links and updates parent link inertia
				transformtoprincipal(childlink);//convert the link's inertia into principal directions first
				idpose.clear();
				assign(childlink,childlink,idpose);//assigns parent link inertia to all the fixed child links.
				assert((childlink->parent_joint->type == Joint::REVOLUTE)||(childlink->parent_joint->type == Joint::PRISMATIC));
				//set name of the joint:
					mbs->joints[index].name = childlink->parent_joint->name;
				//based on type  of joint decide the axis
				switch(childlink->parent_joint->type)
				{
					case Joint::REVOLUTE :
						mbs->joints[index].a.setZero();
						mbs->joints[index].a[0]  = childlink->parent_joint->axis.x;
						mbs->joints[index].a[1]  = childlink->parent_joint->axis.y;
						mbs->joints[index].a[2]  = childlink->parent_joint->axis.z;
						if(childlink->parent_joint->dynamics)
						{
							mbs->damping[index] = childlink->parent_joint->dynamics->damping;
						}
						break;
					case Joint::PRISMATIC :
						mbs->joints[index].a.setZero();
						mbs->joints[index].a[3]  = childlink->parent_joint->axis.x;
						mbs->joints[index].a[4]  = childlink->parent_joint->axis.y;
						mbs->joints[index].a[5]  = childlink->parent_joint->axis.z;
						if(childlink->parent_joint->dynamics)
						{
							mbs->damping[index] = childlink->parent_joint->dynamics->damping;
						}
						break;
				}
				//decide the transformations:
				mbs->joints[index].gp = diffpose(childlink->parent_joint->parent_to_joint_origin_transform,link->inertial->origin);

				//can change to convert to quatxyz2g
				double r,p,y;
				childlink->inertial->origin.rotation.getRPY(r,p,y);

				gcop::SE3::Instance().rpyxyz2g(mbs->joints[index].gc,Vector3d(r,p,y),Vector3d(childlink->inertial->origin.position.x,childlink->inertial->origin.position.y,childlink->inertial->origin.position.z));
				cout<<"level: "<<level<<endl;
				cout<<" Child name "<<childlink->name<<endl;

				index++;
				cout<<"index: "<<index<<endl;

				//create body3d from the current link and make a joint connecting them. This joint is the parent joint of the child link

				//gcop::Body3d<> body;
				//copy m
				mbs->links[index].m = childlink->inertial->mass;
				//copy J
				mbs->links[index].I(0) = childlink->inertial->ixx;
				mbs->links[index].I(1) = childlink->inertial->iyy;
				mbs->links[index].I(2) = childlink->inertial->izz;//we assume the axes are aligned with principal axes
				mbs->links[index].I(3) = childlink->inertial->mass;
				mbs->links[index].I(4) = childlink->inertial->mass;
				mbs->links[index].I(5) = childlink->inertial->mass;
				mbs->links[index].name = childlink->name;
				// make ds = 0;
				mbs->links[index].ds = Vector3d(0,0,0);

				mbs->pis[index] = level; //added parent index to the list.

			}
			walkTree(*child,index,index,mbs);
		}
		else
		{
			for(int j=0;j<level;j++) cout << " "; //indent
			cout << "root link: " << link->name << " has a null child!" << *child << endl;
		}
	}
}
void findnofactivejoints(boost::shared_ptr<Link> link, int &count)
{
	for (vector<boost::shared_ptr<Link> >::const_iterator child = link->child_links.begin(); child != link->child_links.end(); child++)
	{
		if (*child)
		{
			boost::shared_ptr<Link> childlink = *child;
				if((childlink->parent_joint->type == Joint::REVOLUTE)||(childlink->parent_joint->type == Joint::REVOLUTE))
					count++;
			findnofactivejoints(*child,count);
		}
		else
		{
			cout << "root link: " << link->name << " has a null child!" << *child << endl;
		}
			
	}
}
boost::shared_ptr<gcop::Mbs> mbsgenerator(const string &xml_string, string type = "chainbase")
{
	//parse the xml file into a urdf model
	boost::shared_ptr<ModelInterface> urdfmodel = parseURDF(xml_string);

	if (!urdfmodel){
		cerr << "ERROR: Model Parsing the xml failed" << endl;
		return boost::shared_ptr<gcop::Mbs>();
	}
	cout << "robot name is: " << urdfmodel->getName() << endl;

	// get root link
	boost::shared_ptr<Link> root_link = urdfmodel->root_link_;
	if (!root_link) return boost::shared_ptr<gcop::Mbs>();
	//assign 0 inertial to root_link:
	root_link->inertial.reset(new Inertial());


	//create Mbs with the number of links and joints
	int nofactivejoints = 0;
	findnofactivejoints(root_link,nofactivejoints);
	cout<<"Number of active joints: "<<nofactivejoints<<endl;
	boost::shared_ptr<gcop::Mbs> mbs;
	if(type == "airbase")
	{
		mbs.reset(new gcop::Mbs(nofactivejoints+1,4+nofactivejoints));
		mbs->basetype = mbs->AIRBASE;
	}
	else
		mbs.reset(new gcop::Mbs(nofactivejoints+1,6 + nofactivejoints));
	//aggregate and assign rootlink
	root_link->visited = true;//marks the link as visited for aggregation purposes
	Pose idpose;
	idpose.clear();
	aggregate(root_link,root_link,idpose);//aggregates all the fixed joint links and updates parent link inertia
	transformtoprincipal(root_link);
	idpose.clear();
	assign(root_link,root_link,idpose);//assigns parent link inertia to all the fixed child links.

	int index = 0;

	//add the first link to mbs
	mbs->links[index].m = root_link->inertial->mass;
	//copy J
	mbs->links[index].I(0) = root_link->inertial->ixx;
	mbs->links[index].I(1) = root_link->inertial->iyy;
	mbs->links[index].I(2) = root_link->inertial->izz;//we assume the axes are aligned with principal axes
	mbs->links[index].I(3) = root_link->inertial->mass;
	mbs->links[index].I(4) = root_link->inertial->mass;
	mbs->links[index].I(5) = root_link->inertial->mass;
	mbs->links[index].name = root_link->name;
	// make ds = 0;
	mbs->links[index].ds = Vector3d(0,0,0);

	mbs->pis[index] = -1; //added parent index to the list.
	

	cout<<"index: "<<index<<endl;
	cout<<"Level: -1"<<endl;
	cout<<"Root name "<<root_link->name<<endl;
	walkTree(root_link,0,index,mbs);
	mbs->Init();
	//cout<<"Number of links excluding baselink and it's joint: "<<urdfmodel->links_.size()-1 <<"\t"<<urdfmodel->joints_.size()-1<<endl;
	//cout<<"Number of links parsed: "<<es.size()<<endl;
	//cout<<" Number of joints parsed: "<<joints.size()<<endl;
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
    //printf("Could not find the 'robot' element in the xml file");
    model.reset();
    return model;
  }

  // Get robot name
  const char *name = robot_xml->Attribute("name");
  if (!name)
  {
    //printf("No name given for the robot.");
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
        //printf("material '%s' is not unique.", material->name.c_str());
        material.reset();
        model.reset();
        return model;
      }
      else
      {
        model->materials_.insert(make_pair(material->name,material));
        //printf("successfully added a new material '%s'", material->name.c_str());
      }
    }
    else
    {
      //printf("material xml is not initialized correctly");
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
        //printf("link '%s' is not unique.", link->name.c_str());
        model.reset();
        return model;
      }
      else
      {
        // set link visual material
        //printf("setting link '%s' material", link->name.c_str());
        if (link->visual)
        {
          if (!link->visual->material_name.empty())
          {
            if (model->getMaterial(link->visual->material_name))
            {
              //printf("setting link '%s' material to '%s'", link->name.c_str(),link->visual->material_name.c_str());
              link->visual->material = model->getMaterial( link->visual->material_name.c_str() );
            }
            else
            {
              if (link->visual->material)
              {
                //printf("link '%s' material '%s' defined in Visual.", link->name.c_str(),link->visual->material_name.c_str());
                model->materials_.insert(make_pair(link->visual->material->name,link->visual->material));
              }
              else
              {
                //printf("link '%s' material '%s' undefined.", link->name.c_str(),link->visual->material_name.c_str());
                model.reset();
                return model;
              }
            }
          }
        }

        model->links_.insert(make_pair(link->name,link));
        //printf("successfully added a new link '%s'", link->name.c_str());
      }
    }
    else
    {
      //printf("link xml is not initialized correctly");
      model.reset();
      return model;
    }
  }
  if (model->links_.empty()){
    //printf("No link elements found in urdf file");
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
        //printf("joint '%s' is not unique.", joint->name.c_str());
        model.reset();
        return model;
      }
      else
      {
        model->joints_.insert(make_pair(joint->name,joint));
        //printf("successfully added a new joint '%s'", joint->name.c_str());
      }
    }
    else
    {
      //printf("joint xml is not initialized correctly");
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
    //printf("failed to build tree");
    model.reset();
    return model;
  }

  // find the root link
  if (!model->initRoot(parent_link_tree))
  {
    //printf("failed to find root link");
    model.reset();
    return model;
  }
 
  return model;
}

}

