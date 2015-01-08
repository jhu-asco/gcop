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
#include "GLDebugDrawer.h"
#include <stdio.h> //printf debugging
#include <fstream> //Open file

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




VehicleDemo::VehicleDemo()
:
m_cameraHeight(2.7f),
m_minCameraDistance(2.f),
m_maxCameraDistance(10.f),
gVehicleSteering(0),gVehicleVelocity(0),
steeringIncrement(0.04),velocityIncrement(0.1),
tf(5)
{
  m_cameraPosition = btVector3(10,5,10);
  //Create heightField:
  creategreyscaledata();
}

VehicleDemo::~VehicleDemo()
{
  //cleanup in the reverse order of creation/initialization
  delete brccar;
  delete world;
}

void VehicleDemo::initPhysics()
{
  world = new BulletWorld;
  this->m_dynamicsWorld = world->m_dynamicsWorld;

  //btCollisionShape *groundShape = world->CreateGroundPlane(50, 50);
  btCollisionShape *groundShape = world->CreateGroundPlane(20, 20, 0, 20);
  {
    btTransform tr;
    tr.setOrigin(btVector3(0, 0, 0));
    tr.setRotation(btQuaternion(0,0,0));
    world->LocalCreateRigidBody(0,tr, groundShape);
  }

  /////////////Loading Car:
  brccar = new Bulletrccar(*world);//Creating a brand new bullet car system
  brccar->kp_steer = 0.5f;//Steering angle

  //////////// Setup the Optimization Problem:
  cost = new RccarCost(*brccar, gps.Z);
  {
    //Load Params
    Params params;
    params.Load("../../../bin/BulletVehicleEstimation/gnrccarestimation.cfg"); 
    //Gcop Initialization:
    params.GetInt("iters", iters);

    params.GetDouble("tf", tf);
    params.GetInt("N",N);
    assert(tf >0);
    assert(N>0);

    VectorXd R(3);
    if (params.GetVectorXd("R", R))
      cost->R = R.asDiagonal();

    VectorXd S(4);
    if (params.GetVectorXd("S", S))
      cost->S = S.asDiagonal();

    VectorXd P(3);
    if (params.GetVectorXd("P", P)) 
      cost->P = P.asDiagonal();

    cost->UpdateGains();
    //Load Params:
    p0.resize(3);
    if (!params.GetVectorXd("p0", p0)) 
    {
      p0<<2, 50, 0.1;//random values
    }

    params.GetString("sensordata", sensordatafilename);
    params.GetString("controldata", ctrldatafilename);
  }

  //Resize the vectors for states etc

  ////////////Record the sensor data with the true params
  cout<<"Recording Sensor Data"<<endl;
  if(!sensordatafilename.empty() && !ctrldatafilename.empty())
  {
    cout<<"Reading Sensor Data"<<endl;
    std::ifstream ifs1 (ctrldatafilename.c_str(), std::ifstream::in);
    std::ifstream ifs2 (sensordatafilename.c_str(), std::ifstream::in);
    Vector3d z_temp;
    Vector2d us_temp;
    double timez_temp, timeu_temp;
    char comma;
    bool read_result = true;
    N = 0;//Number of segments
    while(ifs1.good())
    {
      if(!(ifs1>>timeu_temp))
        break;
      for(int count1 = 0;count1 < 2;count1++)
      {
        read_result = ifs1.get(comma);
        if(!read_result)
          break;
        read_result = ifs1>>us_temp[count1];
      }

      if(abs(us_temp(0)) < 0.35)//Hacky way of implementing static friction
        us_temp(0) = 0;

      ts.push_back(timeu_temp);
      us.push_back(us_temp);

      cout<<"us["<<(us.size()-1)<<"]: "<<us.back().transpose()<<"\t"<<ts.back()<<endl;
      if(N > 1)
        assert((ts[N] - ts[N-1])>0);
      N++;
    }
    assert(N >=2);//Need atleast two control segments to interpolate the timestep
    ts.push_back(ts[N-1] + (ts[N-1]-ts[N-2]));
    xs.resize(N+1);

    read_result = true;//Reset the state
    int Nz = 0;
    while(ifs2.good())
    {
      if(!(ifs2>>timez_temp))
        break;
      for(int count1 = 0;count1 < 3;count1++)
      {
        read_result = ifs2.get(comma);
        if(!read_result)
          break;
        read_result = ifs2>>z_temp[count1];
      }
      if(!read_result)
        break;
      cout<<timez_temp<<"\t"<<z_temp.transpose()<<endl;
      if(Nz > 0)
      {
        if((timez_temp - ts_sensor[Nz-1]) <= 0)
        {
          cout<<"Skipping Data: New ts[k] is smaller or equal to ts[k-1]; Possible duplicate data: "<<timez_temp<<"\t"<<ts_sensor[Nz-1]<<"\t"<<Nz<<endl;
          continue;
        }
      }
      Nz++;
      ts_sensor.push_back(timez_temp);
      //Convert to Y up and z forward:
      zs.push_back(Vector3d(-z_temp[0], z_temp[2], z_temp[1]));
      cout<<"zs["<<(zs.size()-1)<<"]: "<<zs.back().transpose()<<"\t"<<ts_sensor.back()<<endl;
    }
    cout<<"N: "<<N<<endl;
    cout<<"Nz: "<<Nz<<endl;
    ifs1.close();
    ifs2.close();//Close the streams
    xs[0]<<zs[0][0],zs[0][2],0,0;
    getchar();
  }
  else
  {
    ts.resize(N+1);
    xs.resize(N+1);
    us.resize(N);
    int subdivide = 1;
    zs.resize(N/subdivide);
    ts_sensor.resize(N/subdivide);

    //Set initial posn to 0 and initial velocity to -0.5
    xs[0].setZero();
    //xs[0][3] = -0.5;//testing
    //xs[0][2] = M_PI/4;//testing
    //Initialize the controls and times;
    double h = tf/double(N);

    //Set sensor times:
    for (int k = 0; k <(N/subdivide); ++k)
      ts_sensor[k] = subdivide*k*h;

    ts[0] = 0;
    for(int i = 0;i < N; i++)
    {
      ts[i+1] = (i+1)*h;
      if(i < N/2)
      {
        double temp = (double(2*i)/double(N));
        us[i] = Vector2d(2*temp, 0.2+0.1*temp);
      }
      else
      {
        double temp = (double(2*i - N)/double(N));
        us[i] = Vector2d(2 - 2*temp, 0.3 - 0.3*temp);
      }
    }
    /*for(int i = 0;i < N; i++)
      {
      ts[i+1] = (i+1)*h;
      if(i < N/3)
      {
      double temp = (double(3*i)/double(N));
      us[i] = Vector2d(2*temp, 0.2*temp);
      }
      else if(i < 2*N/3)
      {
      double temp = (double(3*i - N)/double(N));
      us[i] = Vector2d(2 - 1*temp, 0.2 - 0.1*temp);
      }
      else
      {
      double temp = (double(3*i - 2*N)/double(N));
      us[i] = Vector2d(1 - 1 *temp, 0.1 - 0.1*temp);
      }
      }
     */


    clientResetScene();
    int sensor_index = 0;
    //Record the car trajectory
    for (int i =0 ; i < N; i++)
    { 
      if (brccar)
      {
        //cout<<"u: "<<us[i].transpose()<<endl;
        brccar->Step1(xs[i+1], us[i], h);//Move
      }
      if((ts_sensor[sensor_index] - ts[i])>= 0 && (ts_sensor[sensor_index] - ts[i+1]) < 0)
      {
        int near_index = (ts_sensor[sensor_index] - ts[i]) > -(ts_sensor[sensor_index] - ts[i+1])?(i+1):i;
        //Set zsi:
        zs[sensor_index]<<xs[near_index][0],0,xs[near_index][1];
        cout<<"Zs ["<<(sensor_index)<<"]: "<<zs[sensor_index].transpose()<<"ts_sensor: "<<ts_sensor[sensor_index]<<endl;
        sensor_index = sensor_index < (ts_sensor.size()-1)?sensor_index+1:sensor_index;
      }
      //printf("Zs:[%d]: %f\t%f\n",i, zs[i](0), zs[i](2));
    }
  }
  cost->SetReference(&zs, &this->p0);//Set reference for zs

  //Setup the Estimator
  gn = new RccarDoep(*brccar, gps, *cost, ts, xs, us, p0, ts_sensor, &projectmanifold);  


	clientResetScene();
	setCameraDistance(0.5f);
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


void VehicleDemo::MoveAndDisplay(double h)
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

	
	


	
	if (m_dynamicsWorld && brccar)
  {
    Vector2d u;
    u<<gVehicleVelocity, gVehicleSteering;
    //cout<<"u: "<<u.transpose()<<"\tdt: "<<dt<<endl;
    brccar->Step2(u, h);//Move
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
    int size = zs->size();
    for (int i = 0; i < size; i++) {
      const Vector3d &z = (*zs)[i];
      //glVertex3d(-i*0.1, 0.2, -2);
      glVertex3d(z[0], 0.15, z[2]);
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
  gVehicleVelocity = 0.f;
  if(brccar)
  {
    brccar->reset(xs[0]);//Reset the car to initial posn
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
    case 'i' :
      {
        //iterate gauss newton method:
        if(gn)
        {
          gn->debug = false;
          //Display: zs - gn->zs:
          for (int i = 0; i < iters; ++i) {
            //    timer_start(timer);
            gn->Iterate();
            displayCallback();
            //long te = timer_us(timer);
            cout << "Parameter: "<< p0 << endl;
            //cout << "Iteration #" << i << " took: " << te << " us." << endl;
            //getchar();
          }
          for(int i =0;i < N;i++)
          {
            cout<<"Xs["<<i<<"]: "<<xs[i].transpose()<<endl;
          }
          for(int i =0;i < zs.size();i++)
          {
            cout<<"Zs["<<i<<"]: "<<zs[i].transpose()<<endl;
          }
          cout<<p0<<endl;
        }
      } 
    case 'z' :
      {
        clientResetScene();
        for(int i =0 ; i< N;i++)//These use the guessed params
        {
          gVehicleVelocity = us[i](0);
          gVehicleSteering = us[i](1);
          MoveAndDisplay(ts[i+1]-ts[i]);
          usleep((ts[i+1]-ts[i])*1e6);
        }
        displayCallback();
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
