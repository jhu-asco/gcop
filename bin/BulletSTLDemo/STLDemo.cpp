/*
STL Demo adapted from Bullet STL Loader
*/

#include "BulletDynamics/btBulletDynamicsCommon.h"

#include "GLDebugDrawer.h"
#include <stdio.h> //printf debugging

#include "GL_ShapeDrawer.h"

#include "GlutStuff.h"
#include "STLDemo.h"

using namespace gcop;
using namespace std;

const int maxProxies = 32766;
const int maxOverlap = 65535;

STLDemo::STLDemo()
:m_cameraHeight(2.7f),
m_minCameraDistance(0.1f),
m_maxCameraDistance(5.f),
filename(0)
{
  m_cameraPosition = btVector3(10,5,10);
}

STLDemo::~STLDemo()
{
  //Cleanup Bullet class 
  delete world;
}

void STLDemo::initPhysics()
{
  world = new BulletWorld;//Create a Brand new empty Dynamic world
  this->m_dynamicsWorld = world->m_dynamicsWorld;
  if(!filename)
    filename = "../../../bin/BulletSTLDemo/meshes/rocky.stl";
  btCollisionShape *groundShape = world->CreateMeshFromSTL(filename);
  btTransform tr;
  tr.setIdentity();
  tr.setRotation(btQuaternion(0,M_PI/2,0));
  world->LocalCreateRigidBody(0,tr, groundShape);
}


void STLDemo::displayCallback(void) 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
  m_shapeDrawer->drawCoordSystem();

	DemoApplication::renderme();
  

//optional but useful: debug drawing
	if (m_dynamicsWorld)
		m_dynamicsWorld->debugDrawWorld();

	glFlush();
	glutSwapBuffers();
}
