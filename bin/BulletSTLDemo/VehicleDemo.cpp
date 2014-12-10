/*
STL Demo adapted from Bullet STL Loader
*/

#include "BulletDynamics/btBulletDynamicsCommon.h"

#include "GLDebugDrawer.h"
#include <stdio.h> //printf debugging

#include "GL_ShapeDrawer.h"

#include "GlutStuff.h"
#include "VehicleDemo.h"

using namespace gcop;
using namespace std;

const int maxProxies = 32766;
const int maxOverlap = 65535;

VehicleDemo::VehicleDemo()
:m_cameraHeight(2.7f),
m_minCameraDistance(0.1f),
m_maxCameraDistance(5.f),
gVehicleSteering(0),gVehicleVelocity(0),
steeringIncrement(0.04),velocityIncrement(0.1)
{
  m_cameraPosition = btVector3(10,5,10);
}

VehicleDemo::~VehicleDemo()
{
  //Cleanup Bullet class 
  delete world;
}

void VehicleDemo::initPhysics()
{
  world = new BulletWorld;//Create a Brand new empty Dynamic world
  this->m_dynamicsWorld = world->m_dynamicsWorld;
  btCollisionShape *groundShape = world->CreateMeshFromSTL("../../../bin/BulletSTLDemo/meshes/rocky.stl", btVector3(10,10,10));
  btTransform tr;
  tr.setOrigin(btVector3(2, 0, 2));
  tr.setRotation(btQuaternion(0,-M_PI/2,0));
  world->LocalCreateRigidBody(0,tr, groundShape);
  brccar = new Bulletrccar(*world);
}

//to be implemented by the demo
void VehicleDemo::renderme()
{
  updateCamera();

  btScalar m[16];
  int i;


  btVector3 wheelColor(1,0,0);

  btVector3	worldBoundsMin,worldBoundsMax;
  getDynamicsWorld()->getBroadphase()->getBroadphaseAabb(worldBoundsMin,worldBoundsMax);


  if(brccar)
  {
    for (i=0;i<(brccar->m_vehicle)->getNumWheels();i++)
    {
      //synchronize the wheels with the (interpolated) chassis worldtransform
      (brccar->m_vehicle)->updateWheelTransform(i,true);
      //draw wheels (cylinders)
      (brccar->m_vehicle)->getWheelInfo(i).m_worldTransform.getOpenGLMatrix(m);
      m_shapeDrawer->drawOpenGL(m,(brccar->m_wheelShape),wheelColor,getDebugMode(),worldBoundsMin,worldBoundsMax);
    }
  }

  DemoApplication::renderme();

}

void VehicleDemo::clientMoveAndDisplay(void)
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	
	if (m_dynamicsWorld && brccar)
  {
    Vector2d u;
    u<<gVehicleVelocity, gVehicleSteering;
    float dt = getDeltaTimeMicroseconds() * 0.000001f;

    //cout<<"u: "<<u.transpose()<<"\tdt: "<<dt<<endl;
    brccar->Step2(u, dt);//Move //can change number of substeps based on idle #TODO
  }


#ifdef USE_QUICKPROF 
        btProfiler::beginBlock("render"); 
#endif //USE_QUICKPROF 


	renderme(); 

	//optional but useful: debug drawing
	if (m_dynamicsWorld)
		m_dynamicsWorld->debugDrawWorld();

#ifdef USE_QUICKPROF 
        btProfiler::endBlock("render"); 
#endif 
	

	glFlush();
	glutSwapBuffers();

}

void VehicleDemo::displayCallback(void) 
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

void VehicleDemo::clientResetScene()
{
  gVehicleSteering = 0.f;
  gVehicleVelocity = 0.f;
  if(brccar)
  {
    Vector4d x0;
    brccar->initialz = 2;
    x0<<0,0,M_PI,0;
    brccar->reset(x0);//Reset the car to initial posn
		for (int i=0;i<(brccar->m_vehicle->getNumWheels());i++)
		{
			//synchronize the wheels with the (interpolated) chassis worldtransform
			(brccar->m_vehicle)->updateWheelTransform(i,true);
		}
	}
}

void VehicleDemo::specialKeyboardUp(int key, int x, int y)
{
//	printf("keyup = %i x=%i y=%i\n",key,x,y);
   switch (key) 
    {
    case GLUT_KEY_UP :
		{
			//gVehicleVelocity = 0.f;
		break;
		}
	case GLUT_KEY_DOWN :
		{			
			//gVehicleVelocity = 0.f; 
		break;
		}
	default:
		DemoApplication::specialKeyboardUp(key,x,y);
        break;
    }

}

void VehicleDemo::keyboardCallback(unsigned char key, int x, int y)
{
  printf("keyup = %i x=%i y=%i\n",key,x,y);

  switch(key)
  {
    case 'r' :
      {
        clientResetScene();
        break;
      }
    default:
      {
        DemoApplication::keyboardCallback(key, x, y);
      }
  }
}

void VehicleDemo::specialKeyboard(int key, int x, int y)
{

//	printf("key = %i x=%i y=%i\n",key,x,y);

    switch (key) 
    {
      case GLUT_KEY_LEFT : 
        {
          if(brccar)
          {
          gVehicleSteering += steeringIncrement;
            if (	gVehicleSteering > brccar->steeringClamp)
              gVehicleSteering = brccar->steeringClamp;
          }

          break;
        }
      case GLUT_KEY_RIGHT : 
        {
          if(brccar)
          {
            gVehicleSteering -= steeringIncrement;
            if (	gVehicleSteering < -(brccar->steeringClamp))
              gVehicleSteering = -(brccar->steeringClamp);
          }

          break;
        }
      case GLUT_KEY_UP :
        {
          if(brccar)
          {
            gVehicleVelocity += velocityIncrement;
            if( gVehicleVelocity > (brccar->velocityClamp))
              gVehicleVelocity = (brccar->velocityClamp);
          }

          break;
        }
      case GLUT_KEY_DOWN :
        {			
          if(brccar)
          {
            gVehicleVelocity -= velocityIncrement;
            if( gVehicleVelocity < - (brccar->velocityClamp))
              gVehicleVelocity = -(brccar->velocityClamp);
          }

          break;
        } 
      default:
        DemoApplication::specialKeyboard(key,x,y);
        break;
    }

//	glutPostRedisplay();
}

void	VehicleDemo::updateCamera()
{
	
//#define DISABLE_CAMERA 1
#ifdef DISABLE_CAMERA
	DemoApplication::updateCamera();
	return;
#endif //DISABLE_CAMERA

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	btTransform chassisWorldTrans;

	//look at the vehicle
  if(brccar)
  {
    chassisWorldTrans = brccar->m_carChassis->getWorldTransform();
    m_cameraTargetPosition = chassisWorldTrans.getOrigin();
  }

	//interpolate the camera height
#ifdef FORCE_ZAXIS_UP
	m_cameraPosition[2] = (15.0*m_cameraPosition[2] + m_cameraTargetPosition[2] + m_cameraHeight)/16.0;
#else
	m_cameraPosition[1] = (15.0*m_cameraPosition[1] + m_cameraTargetPosition[1] + m_cameraHeight)/16.0;
#endif 

	btVector3 camToObject = m_cameraTargetPosition - m_cameraPosition;

	//keep distance between min and max distance
	float cameraDistance = camToObject.length();
	float correctionFactor = 0.f;
	if (cameraDistance < m_minCameraDistance)
	{
		correctionFactor = 0.15*(m_minCameraDistance-cameraDistance)/cameraDistance;
	}
	if (cameraDistance > m_maxCameraDistance)
	{
		correctionFactor = 0.15*(m_maxCameraDistance-cameraDistance)/cameraDistance;
	}
	m_cameraPosition -= correctionFactor*camToObject;
	
	  btScalar aspect = m_glutScreenWidth / (btScalar)m_glutScreenHeight;
        glFrustum (-aspect, aspect, -1.0, 1.0, 1.0, 10000.0);

         glMatrixMode(GL_MODELVIEW);
         glLoadIdentity();

    gluLookAt(m_cameraPosition[0],m_cameraPosition[1],m_cameraPosition[2],
                      m_cameraTargetPosition[0],m_cameraTargetPosition[1], m_cameraTargetPosition[2],
                          m_cameraUp.getX(),m_cameraUp.getY(),m_cameraUp.getZ());
}
