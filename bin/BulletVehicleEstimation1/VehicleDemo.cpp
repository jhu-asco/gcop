/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  http://continuousphysics.com/Bullet/

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

/// September 2006: VehicleDemo is work in progress, this file is mostly just a placeholder
/// This VehicleDemo file is very early in development, please check it later
/// One todo is a basic engine model:
/// A function that maps user input (throttle) into torque/force applied on the wheels
/// with gears etc.
#include "BulletDynamics/btBulletDynamicsCommon.h"
#include "BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h"
#include "helper_header.h"
//extern char MyHeightfield[];
//
// By default, Bullet Vehicle uses Y as up axis.
// You can override the up axis, for example Z-axis up. Enable this define to see how to:
//#define FORCE_ZAXIS_UP 1
//

#ifdef FORCE_ZAXIS_UP
		int rightIndex = 0; 
		int upIndex = 2; 
		int forwardIndex = 1;
		btVector3 wheelDirectionCS0(0,0,-1);
		btVector3 wheelAxleCS(1,0,0);
#else
		int rightIndex = 0;
		int upIndex = 1;
		int forwardIndex = 2;
		btVector3 wheelDirectionCS0(0,-1,0);
		btVector3 wheelAxleCS(-1,0,0);
#endif

#include "GLDebugDrawer.h"
#include <stdio.h> //printf debugging

#include "GL_ShapeDrawer.h"

#include "GlutStuff.h"
#include "VehicleDemo.h"

using namespace gcop;
using namespace Eigen;
using namespace std;

const int maxProxies = 32766;
const int maxOverlap = 65535;

///btRaycastVehicle is the interface for the constraint that implements the raycast vehicle
///notice that for higher-quality slow-moving vehicles, another approach might be better
///implementing explicit hinged-wheel constraints with cylinder collision, rather then raycasts
float	gEngineForce = 0.f;
float	gBreakingForce = 0.f;

float	maxEngineForce = 25.f;//this should be engine/velocity dependent
float	maxBreakingForce = 2.f;

float	gVehicleSteering = 0.f;
float	steeringIncrement = 0.04f;
float	steeringClamp = 0.3f;
float	wheelRadius = 0.07f;
float	wheelWidth = 0.06f;
//float	wheelFriction = 200;//BT_LARGE_FLOAT;
float	wheelFriction = BT_LARGE_FLOAT;
float	suspensionStiffness = 1000.f;
float	suspensionDamping = 100.f;
float	suspensionCompression = 0.02f;
float	rollInfluence = 0.01f;//1.0f;


btScalar suspensionRestLength(0.1);

#define CAR_HALFLENGTH 0.32
#define CAR_HALFWIDTH 0.19
#define CAR_HALFHEIGHT 0.1



////////////////////////////////////




VehicleDemo::VehicleDemo()
:
m_carChassis(0),
m_indexVertexArrays(0),
m_vertices(0),
m_cameraHeight(2.7f),
m_minCameraDistance(0.1f),
//m_maxCameraDistance(5.f),
m_maxCameraDistance(10.f),
tf(5)
//sys(), gps(),
//cost(sys,gps.Z)
{
  m_vehicle = 0;
  m_wheelShape = 0;
  m_cameraPosition = btVector3(10,5,10);
  //Create heightField:
  creategreyscaledata();
  cout<<"1"<<endl;
  cost = new RccarCost(sys, gps.Z);
  cout<<"2"<<endl;
  //Load Params
  Params params;
  params.Load("../../../bin/BulletVehicleEstimation1/gnrccarestimation.cfg"); 
  //Gcop Initialization:
  params.GetInt("iters", iters);

  params.GetDouble("tf", tf);
  params.GetInt("N",N);
  assert(tf >0);
  assert(N>0);

  {
    VectorXd R(3);
    if (params.GetVectorXd("R", R))
      cost->R = R.asDiagonal();

    VectorXd S(4);
    if (params.GetVectorXd("S", S))
      cost->S = S.asDiagonal();

    VectorXd P(2);
    if (params.GetVectorXd("P", P)) 
      cost->P = P.asDiagonal();

    cost->UpdateGains();
  }

  //Resize the vectors for states etc
  ts.resize(N+1);
  xs.resize(N+1);
  zs.resize(N);
  us.resize(N);
  p0.resize(2);
  //Initialize the controls and times;
  double h = tf/double(N);
  ts[0] = 0;
  for(int i = 0;i < N; i++)
  {
    ts[i+1] = (i+1)*h;
    if(i < N/3)
      us[i] = Vector2d(2, 0.2*(double(3*i)/double(N)));
    else if(i < 2*N/3)
      us[i] = Vector2d(3, 0.2 - 0.3*(double(3*i - N)/double(N)));//Verify if steering inputs are torques or angles?
      //us[i] = Vector2d(-8, -0.2);//Verify if steering inputs are torques or angles?
    else
      us[i] = Vector2d(-4, -0.1 - 0.2*(double(3*i - 2*N)/double(N)));
  }
  //Load Params:
  if (!params.GetVectorXd("p0", p0)) 
  {
    p0<<0.5, 0.01;//random values
  }
  gn = new RccarDoep(sys, gps, *cost, ts, xs, us, p0, &projectmanifold);  
}

VehicleDemo::~VehicleDemo()
{
		//cleanup in the reverse order of creation/initialization

	//remove the rigidbodies from the dynamics world and delete them
	int i;
	for (i=m_dynamicsWorld->getNumCollisionObjects()-1; i>=0 ;i--)
	{
		btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if (body && body->getMotionState())
		{
			delete body->getMotionState();
		}
		m_dynamicsWorld->removeCollisionObject( obj );
		delete obj;
	}

	//delete collision shapes
	for (int j=0;j<m_collisionShapes.size();j++)
	{
		btCollisionShape* shape = m_collisionShapes[j];
		delete shape;
	}

	delete m_indexVertexArrays;
	delete m_vertices;

	//delete dynamics world
	delete m_dynamicsWorld;

	delete m_vehicleRayCaster;

	delete m_vehicle;

	delete m_wheelShape;

	//delete solver
	delete m_constraintSolver;

	//delete broadphase
	delete m_overlappingPairCache;

	//delete dispatcher
	delete m_dispatcher;

	delete m_collisionConfiguration;

}

void VehicleDemo::initPhysics()
{
	
#ifdef FORCE_ZAXIS_UP
	m_cameraUp = btVector3(0,0,1);
	m_forwardAxis = 1;
#endif

	btCollisionShape* groundShape = new btBoxShape(btVector3(50,3,50));
	m_collisionShapes.push_back(groundShape);
	m_collisionConfiguration = new btDefaultCollisionConfiguration();
	m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);
	btVector3 worldMin(-1000,-1000,-1000);
	btVector3 worldMax(1000,1000,1000);
	m_overlappingPairCache = new btAxisSweep3(worldMin,worldMax);
	m_constraintSolver = new btSequentialImpulseConstraintSolver();
	m_dynamicsWorld = new btDiscreteDynamicsWorld(m_dispatcher,m_overlappingPairCache,m_constraintSolver,m_collisionConfiguration);
#ifdef FORCE_ZAXIS_UP
	m_dynamicsWorld->setGravity(btVector3(0,0,-10));
#endif 

	//m_dynamicsWorld->setGravity(btVector3(0,0,0));
btTransform tr;
tr.setIdentity();

//either use heightfield or triangle mesh
#define  USE_TRIMESH_GROUND 1
#ifdef USE_TRIMESH_GROUND
	int i;

const float TRIANGLE_SIZE=20.f;

	//create a triangle-mesh ground
	int vertStride = sizeof(btVector3);
	int indexStride = 3*sizeof(int);

	const int NUM_VERTS_X = 20;
	const int NUM_VERTS_Y = 20;
	const int totalVerts = NUM_VERTS_X*NUM_VERTS_Y;
	
	const int totalTriangles = 2*(NUM_VERTS_X-1)*(NUM_VERTS_Y-1);

	m_vertices = new btVector3[totalVerts];
	int*	gIndices = new int[totalTriangles*3];

	

	for ( i=0;i<NUM_VERTS_X;i++)
	{
		for (int j=0;j<NUM_VERTS_Y;j++)
		{
			float wl = .2f;
			//height set to zero, but can also use curved landscape, just uncomment out the code
			float height = 0.f;//20.f*sinf(float(i)*wl)*cosf(float(j)*wl);
#ifdef FORCE_ZAXIS_UP
			m_vertices[i+j*NUM_VERTS_X].setValue(
				(i-NUM_VERTS_X*0.5f)*TRIANGLE_SIZE,
				(j-NUM_VERTS_Y*0.5f)*TRIANGLE_SIZE,
				height
				);

#else
			m_vertices[i+j*NUM_VERTS_X].setValue(
				(i-NUM_VERTS_X*0.5f)*TRIANGLE_SIZE,
				height,
				(j-NUM_VERTS_Y*0.5f)*TRIANGLE_SIZE);
#endif

		}
	}

	int index=0;
	for ( i=0;i<NUM_VERTS_X-1;i++)
	{
		for (int j=0;j<NUM_VERTS_Y-1;j++)
		{
			gIndices[index++] = j*NUM_VERTS_X+i;
			gIndices[index++] = j*NUM_VERTS_X+i+1;
			gIndices[index++] = (j+1)*NUM_VERTS_X+i+1;

			gIndices[index++] = j*NUM_VERTS_X+i;
			gIndices[index++] = (j+1)*NUM_VERTS_X+i+1;
			gIndices[index++] = (j+1)*NUM_VERTS_X+i;
		}
	}
	
	m_indexVertexArrays = new btTriangleIndexVertexArray(totalTriangles,
		gIndices,
		indexStride,
		totalVerts,(btScalar*) &m_vertices[0].x(),vertStride);

	bool useQuantizedAabbCompression = true;
	groundShape = new btBvhTriangleMeshShape(m_indexVertexArrays,useQuantizedAabbCompression);
	
	//tr.setOrigin(btVector3(0,-4.5f,0));
	tr.setOrigin(btVector3(0,-0.15f,0));

#else
	//testing btHeightfieldTerrainShape
	//int width=128;
	//int height=128;
	

#ifdef LOAD_FROM_FILE
	unsigned char* heightfieldData = new unsigned char[width*height];
	{
		for (int i=0;i<width*height;i++)
		{
			heightfieldData[i]=0;
		}
	}

	char*	filename="heightfield128x128.raw";
	FILE* heightfieldFile = fopen(filename,"r");
	if (!heightfieldFile)
	{
		filename="../../heightfield128x128.raw";
		heightfieldFile = fopen(filename,"r");
	}
	if (heightfieldFile)
	{
		int numBytes =fread(heightfieldData,1,width*height,heightfieldFile);
		//btAssert(numBytes);
		if (!numBytes)
		{
			printf("couldn't read heightfield at %s\n",filename);
		}
		fclose (heightfieldFile);
	}
#else
	char* heightfieldData = MyHeightfield;
#endif


	//btScalar maxHeight = 20000.f;//exposes a bug
	//btScalar maxHeight = 100;
	btScalar maxHeight = 128;
	
	bool useFloatDatam=false;
	bool flipQuadEdges=false;

	btHeightfieldTerrainShape* heightFieldShape = new btHeightfieldTerrainShape(width,height,heightfieldData,maxHeight,upIndex,useFloatDatam,flipQuadEdges);
	btVector3 mmin,mmax;
	heightFieldShape->getAabb(btTransform::getIdentity(),mmin,mmax);

	groundShape = heightFieldShape;
	
	heightFieldShape->setUseDiamondSubdivision(true);

	//btVector3 localScaling(100,1,100);
	btVector3 localScaling(0.1,0.1,0.1);
	localScaling[upIndex]=0.5f;
	groundShape->setLocalScaling(localScaling);

	//tr.setOrigin(btVector3(0,9940,0));
	//tr.setOrigin(btVector3(0,49.4,0));
  tr.setOrigin(btVector3(0,30,0));
  printf("hello\n");

#endif //

	m_collisionShapes.push_back(groundShape);

	//create ground object
	localCreateRigidBody(0,tr,groundShape);
	tr.setOrigin(btVector3(0,0,0));//-64.5f,0));

#ifdef FORCE_ZAXIS_UP
//   indexRightAxis = 0; 
//   indexUpAxis = 2; 
//   indexForwardAxis = 1; 
	btCollisionShape* chassisShape = new btBoxShape(btVector3(1.f,2.f, 0.5f));
	btCompoundShape* compound = new btCompoundShape();
	btTransform localTrans;
	 //may be on the right track now to fixing that. 
   //Heavier masses and normal grav make the steering much more responsive with less tendency to move thelocalTrans.setIdentity();
	//localTrans effectively shifts the center of mass with respect to the chassis
	localTrans.setOrigin(btVector3(0,0,1));
#else
	btCollisionShape* chassisShape = new btBoxShape(btVector3(CAR_HALFWIDTH,CAR_HALFHEIGHT,CAR_HALFLENGTH));
	m_collisionShapes.push_back(chassisShape);

	btCompoundShape* compound = new btCompoundShape();
	m_collisionShapes.push_back(compound);
	btTransform localTrans;
	localTrans.setIdentity();
	//localTrans effectively shifts the center of mass with respect to the chassis
  localTrans.setOrigin(btVector3(0,0.0,0));
	//localTrans.setOrigin(btVector3(0,0.5,0));
#endif

	compound->addChildShape(localTrans,chassisShape);

	tr.setOrigin(btVector3(0,0.0f,0));

	//m_carChassis = localCreateRigidBody(800,tr,compound);//chassisShape);
	m_carChassis = localCreateRigidBody(20,tr,compound);//chassisShape);
	//m_carChassis->setDamping(0.2,0.2);
	
	m_wheelShape = new btCylinderShapeX(btVector3(wheelWidth,wheelRadius,wheelRadius));
	
	clientResetScene();

	/// create vehicle
	{
		
		m_vehicleRayCaster = new btDefaultVehicleRaycaster(m_dynamicsWorld);
		m_vehicle = new btRaycastVehicle(m_tuning,m_carChassis,m_vehicleRayCaster);
		
		///never deactivate the vehicle
		m_carChassis->setActivationState(DISABLE_DEACTIVATION);

		m_dynamicsWorld->addVehicle(m_vehicle);

		float connectionHeight = 0.0f + 0.17f - 0.03f - CAR_HALFHEIGHT;// Localtrans + Suspension tip connection height - ground clearance - half the height of car

	
		bool isFrontWheel=true;

		//choose coordinate system
		m_vehicle->setCoordinateSystem(rightIndex,upIndex,forwardIndex);

#ifdef FORCE_ZAXIS_UP
		btVector3 connectionPointCS0(CUBE_HALF_EXTENTS-(0.3*wheelWidth),2*CUBE_HALF_EXTENTS-wheelRadius, connectionHeight);
#else
		btVector3 connectionPointCS0(0.19,connectionHeight,CAR_HALFLENGTH-wheelRadius);
#endif

		m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);
#ifdef FORCE_ZAXIS_UP
		connectionPointCS0 = btVector3(-CUBE_HALF_EXTENTS+(0.3*wheelWidth),2*CUBE_HALF_EXTENTS-wheelRadius, connectionHeight);
#else
		connectionPointCS0 = btVector3(-0.19,connectionHeight,CAR_HALFLENGTH-wheelRadius);
#endif

		m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);
#ifdef FORCE_ZAXIS_UP
		connectionPointCS0 = btVector3(-CUBE_HALF_EXTENTS+(0.3*wheelWidth),-2*CUBE_HALF_EXTENTS+wheelRadius, connectionHeight);
#else
		connectionPointCS0 = btVector3(-0.19,connectionHeight,-CAR_HALFLENGTH+wheelRadius);
#endif //FORCE_ZAXIS_UP
		isFrontWheel = false;
		m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);
#ifdef FORCE_ZAXIS_UP
		connectionPointCS0 = btVector3(CUBE_HALF_EXTENTS-(0.3*wheelWidth),-2*CUBE_HALF_EXTENTS+wheelRadius, connectionHeight);
#else
		connectionPointCS0 = btVector3(0.19,connectionHeight,-CAR_HALFLENGTH+wheelRadius);
#endif
		m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);
		
		for (int i=0;i<m_vehicle->getNumWheels();i++)
		{
			btWheelInfo& wheel = m_vehicle->getWheelInfo(i);
			wheel.m_suspensionStiffness = suspensionStiffness;
			wheel.m_wheelsDampingRelaxation = suspensionDamping;
			wheel.m_wheelsDampingCompression = suspensionCompression;
			wheel.m_frictionSlip = wheelFriction;
			wheel.m_rollInfluence = rollInfluence;
		}
	}

	
	setCameraDistance(0.5f);

}


void VehicleDemo::sampleSensor(Vector3d &z)
{
  if(m_carChassis && m_vehicle)
  {
    btVector3 &vec = (m_carChassis->getWorldTransform()).getOrigin();
    z[0] = -vec[2];
    z[1] = -vec[0];//NOTE: Coord systems different
    z[2] = 0;//Planar sampling for now
  }
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



	for (i=0;i<m_vehicle->getNumWheels();i++)
	{
		//synchronize the wheels with the (interpolated) chassis worldtransform
		m_vehicle->updateWheelTransform(i,true);
		//draw wheels (cylinders)
		m_vehicle->getWheelInfo(i).m_worldTransform.getOpenGLMatrix(m);
		m_shapeDrawer->drawOpenGL(m,m_wheelShape,wheelColor,getDebugMode(),worldBoundsMin,worldBoundsMax);
	}


	DemoApplication::renderme();

}

void VehicleDemo::clientMoveAndDisplay()
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

	
	{			
		int wheelIndex = 2;
		m_vehicle->applyEngineForce(gEngineForce,wheelIndex);
		m_vehicle->setBrake(gBreakingForce,wheelIndex);
		wheelIndex = 3;
		m_vehicle->applyEngineForce(gEngineForce,wheelIndex);
		m_vehicle->setBrake(gBreakingForce,wheelIndex);


		wheelIndex = 0;
		m_vehicle->setSteeringValue(gVehicleSteering,wheelIndex);
		wheelIndex = 1;
		m_vehicle->setSteeringValue(gVehicleSteering,wheelIndex);
	}


	float dt = getDeltaTimeMicroseconds() * 0.000001f;
	
	if (m_dynamicsWorld)
	{
		//during idle mode, just run 1 simulation step maximum
		if (m_idle)
			dt = tf/double(N);
		int maxSimSubSteps = m_idle ? round(500*dt) : 2;

		int numSimSteps = m_dynamicsWorld->stepSimulation(dt,maxSimSubSteps, 1.0f/500.0f);

//#define VERBOSE_FEEDBACK
#ifdef VERBOSE_FEEDBACK
		if (!numSimSteps)
			printf("Interpolated transforms\n");
		else
		{
			if (numSimSteps > maxSimSubSteps)
			{
				//detect dropping frames
				printf("Dropped (%i) simulation steps out of %i\n",numSimSteps - maxSimSubSteps,numSimSteps);
			} else
			{
				printf("Simulated (%i) steps\n",numSimSteps);
			}
		}
#endif //VERBOSE_FEEDBACK

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

void VehicleDemo::renderTrajectory(vector<Vector3d> *zs, btVector3 *color)
{
  if(zs)
  {
    glLineWidth(5);
    if(color)
      glColor3f(color->x(), color->y(), color->z());
    else
      glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < N; i++) {
      const Vector3d &z = (*zs)[i];
      //glVertex3d(-i*0.1, 0.2, -2);
      glVertex3d(-z[1], 0.2, -z[0]);
    }
    glEnd();
  }
}



void VehicleDemo::displayCallback(void) 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
  m_shapeDrawer->drawCoordSystem();

  //resetPerspectiveProjection();
  //glPopMatrix();
  btVector3 color(1,0,0);
  renderTrajectory(&(this->zs), &color);
  color = btVector3(0,1,0);
  renderTrajectory(&(gn->zs), &color);

	renderme();
  

//optional but useful: debug drawing
	if (m_dynamicsWorld)
		m_dynamicsWorld->debugDrawWorld();

	glFlush();
	glutSwapBuffers();
}



void VehicleDemo::clientResetScene()
{
	gVehicleSteering = 0.f;
	m_carChassis->setCenterOfMassTransform(btTransform::getIdentity());
	m_carChassis->setLinearVelocity(btVector3(0,0,0));
	m_carChassis->setAngularVelocity(btVector3(0,0,0));
	m_dynamicsWorld->getBroadphase()->getOverlappingPairCache()->cleanProxyFromPairs(m_carChassis->getBroadphaseHandle(),getDynamicsWorld()->getDispatcher());
	if (m_vehicle)
	{
		m_vehicle->resetSuspension();
		for (int i=0;i<m_vehicle->getNumWheels();i++)
		{
			//synchronize the wheels with the (interpolated) chassis worldtransform
			m_vehicle->updateWheelTransform(i,true);
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
			gEngineForce = 0.f;
		break;
		}
	case GLUT_KEY_DOWN :
		{			
			gBreakingForce = 0.f; 
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
    case 'i' :
      {
        //iterate gauss newton method:
        if(gn)
        {
          gn->debug = false;

          for (int i = 0; i < iters; ++i) {
            //    timer_start(timer);
            gn->Iterate();
            displayCallback();
            //long te = timer_us(timer);
            cout << "Parameter: "<< p0 << endl;
            //cout << "Iteration #" << i << " took: " << te << " us." << endl;
            getchar();
          }
        }
      }
    case 'z' :
      {
        this->clientResetScene();//Reset the scene
        //Record the car trajectory
        for (int i =0 ; i < N; i++)
        { 
          gEngineForce = -us[i](0);
          gBreakingForce = 0.f;
          gVehicleSteering = -us[i](1);
          m_idle = true;
          this->clientMoveAndDisplay();
          //btMatrix3x3 wheelbasis = (m_carChassis->getWorldTransform().inverse()* m_vehicle->getWheelTransformWS(0)).getBasis();
          //btVector3 wheelrpybasis;
          //wheelbasis.getEulerZYX(wheelrpybasis[2],wheelrpybasis[1],wheelrpybasis[0]);
          this->sampleSensor(zs[i]);
          //printf("Zs:[%d]: %f\t%f\t Wheel_angle: %f\t%f\t%f\n",i, zs[i](0), zs[i](1), wheelrpybasis[0], wheelrpybasis[1], wheelrpybasis[2]);
          printf("Zs:[%d]: %f\t%f\n",i, zs[i](0), zs[i](1));
        }
        //set this to the initial state:
        xs[0].setZero();
        xs[0].head(2) = zs[0].head(2);
        //Update xs, us and zs:
        gn->Update(false);
        //Initial Trajectory:
        for(int i=0; i<=N;i++)
        {
          std::cout<<"Xs["<<i<<"]: "<<xs[i].transpose()<<std::endl;
        }
        cost->SetReference(&zs, &this->p0);//Set reference for zs
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
          gVehicleSteering += steeringIncrement;
          if (	gVehicleSteering > steeringClamp)
            gVehicleSteering = steeringClamp;

          break;
        }
      case GLUT_KEY_RIGHT : 
        {
          gVehicleSteering -= steeringIncrement;
          if (	gVehicleSteering < -steeringClamp)
            gVehicleSteering = -steeringClamp;

          break;
        }
      case GLUT_KEY_UP :
        {
          gEngineForce = maxEngineForce;
          gBreakingForce = 0.f;
          break;
        }
      case GLUT_KEY_DOWN :
        {			
          gBreakingForce = maxBreakingForce; 
          gEngineForce = 0.f;
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
	m_carChassis->getMotionState()->getWorldTransform(chassisWorldTrans);
	m_cameraTargetPosition = chassisWorldTrans.getOrigin();

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
