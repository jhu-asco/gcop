#include "bulletrccar.h"

using namespace gcop;
using namespace Eigen;

Bulletrccar::Bulletrccar(btDynamicsWorld*		m_dynamicsWorld_,btAlignedObjectArray<btCollisionShape*> &m_collisionShapes, DemoApplication* appinstance) : 
  Rccar(3),m_carChassis(0), m_vehicle(0), m_wheelShape(0)
  ,carmass(20.0f),rightIndex(0),upIndex(1),forwardIndex(2)
  ,gEngineForce(0),gBreakingForce(0.f),maxEngineForce(100.f),maxBreakingForce(10.f)
  ,gVehicleSteering(0.f),steeringClamp(0.5f), velocityClamp(5.f)
  ,wheelRadius(0.07f), wheelWidth(0.06f),wheelFriction(BT_LARGE_FLOAT)
  ,suspensionStiffness(1000.f),suspensionDamping(100.f),suspensionCompression(0.02f)
  ,rollInfluence(0.1f),suspensionRestLength(0.122)
  ,m_defaultContactProcessingThreshold(BT_LARGE_FLOAT)
  , m_dynamicsWorld(m_dynamicsWorld_)
  , gain_cmdvelocity(2), kp_torque(50), kp_steer(0.1)
{
  wheelDirectionCS0.setValue(0,-1,0);
  wheelAxleCS.setValue(-1,0,0);
  this->car_halfdims[0] = 0.19f;
  this->car_halfdims[1] = 0.1f;
  this->car_halfdims[2] =  0.32f;
  {
    //The car details are hardcoded
    btCollisionShape* chassisShape = new btBoxShape(btVector3(car_halfdims[0],car_halfdims[1],car_halfdims[2]));//Need to use half extents of the box to create the chassis shape  in BULLET !!	
    m_collisionShapes.push_back(chassisShape);

    {
      btTransform initialtr;//temporary transform
      initialtr.setIdentity();
      m_carChassis = appinstance->localCreateRigidBody(carmass,initialtr,chassisShape);//(chassis Rigid Body)
    }

    //Loading Wheel Shapes:
    m_wheelShape = new btCylinderShapeX(btVector3(wheelWidth,wheelRadius,wheelRadius));

    /// create vehicle
    {
      // Ray casting helper class for the vehicle
      m_vehicleRayCaster = new btDefaultVehicleRaycaster(m_dynamicsWorld);
      // Ray Casting Vehicle is made from existing rigid body. It converts a rigid body into a vehicle by attaching four rays looking downwards at four ends which are equivalent to wheels and uses the forces generated to move the rigid body.
      m_vehicle = new btRaycastVehicle(m_tuning,m_carChassis,m_vehicleRayCaster);

      ///never deactivate the vehicle from computing collisions etc. Usually the bodies are deactivated after some time if they are not moving
      m_carChassis->setActivationState(DISABLE_DEACTIVATION);

      m_dynamicsWorld->addVehicle(m_vehicle);

      float connectionHeight = 0.17f - 0.03f - car_halfdims[1];// Suspension tip connection height - ground clearance - half the height of car


      bool isFrontWheel=true;
      //choose coordinate system
      m_vehicle->setCoordinateSystem(rightIndex,upIndex,forwardIndex);

      //Front Right Wheel
      btVector3 connectionPointCS0(0.19,connectionHeight,car_halfdims[2]-wheelRadius);

      m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

      //Front Left Wheel
      connectionPointCS0 = btVector3(-0.19,connectionHeight,car_halfdims[2]-wheelRadius);
      m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

      //Back Right Wheel
      connectionPointCS0 = btVector3(-0.19,connectionHeight,-car_halfdims[2]+wheelRadius);

      isFrontWheel = false;
      m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

      //Back Left Wheel
      connectionPointCS0 = btVector3(0.19,connectionHeight,-car_halfdims[2]+wheelRadius);
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
  }

  U.bnd = true;
  U.lb[0] = -velocityClamp;
  U.ub[0] = velocityClamp;//Max cmd
  U.lb[1] = -steeringClamp;
  U.ub[1] = steeringClamp;
}


double Bulletrccar::Step1(Vector4d& xb, const Vector2d& u, 
                   double h, const VectorXd *p,
                   Matrix4d *A, Matrix42d *B, Matrix4pd *C) {
  //Set Parameters
  if(p)
  {
    gain_cmdvelocity = (*p)[0];

    if((*p)[1] > 0)
      kp_torque = (*p)[1];

    if((*p)[2] < 0)  
      kp_steer = 0;
    else if((*p)[2] > 1)  
      kp_steer = 1;
    else
      kp_steer = (*p)[2];
  }
  //Applies the Forces to the Vehicle
  {
    //Set Controls
    gEngineForce = kp_torque*(gain_cmdvelocity*u[0] - this->x(3));
    gVehicleSteering = (1-kp_steer)*gVehicleSteering + kp_steer*u[1];
    //std::cout<<"Engine Force: "<<gEngineForce<<"\tgVehicleSteering: "<<gVehicleSteering<<std::endl;
    if(!m_vehicle)
    {
      cout<<"M_Vehicle not defined"<<endl;
    }
    /*if(abs(engineForce) > 0 && abs(breakingForce) > 0)
      {
      cerr<<"Both Engine and Breaking forces are nonzero"<<endl;
      return;
      }
      m_vehicle->setBrake(breakingForce,wheelIndex);
      m_vehicle->setBrake(breakingForce,wheelIndex);
     */
    int wheelIndex = 2;
    m_vehicle->applyEngineForce(gEngineForce,wheelIndex);

    wheelIndex = 3;
    m_vehicle->applyEngineForce(gEngineForce,wheelIndex);
    wheelIndex = 0;
    m_vehicle->setSteeringValue(gVehicleSteering,wheelIndex);
    wheelIndex = 1;
    m_vehicle->setSteeringValue(gVehicleSteering,wheelIndex);
  }

  //Set Linearization parameters:
  const double v = this->x[3];
  double c = cos(this->x[2]);
  double s = sin(this->x[2]);//Can find a better derivative based on initial and final values #TODO
  if (A) {
    A->setIdentity();
    (*A)(0,2) = -h*s*v;
    (*A)(0,3) = h*c;
    (*A)(1,2) = h*c*v;
    (*A)(1,3) = h*s;
    (*A)(2,3) = h*tan(gVehicleSteering)/l;
  }

  if (B) {
    B->setZero();
    (*B)(2,1) =  h*kp_steer/(l*v*pow(cos(gVehicleSteering),2));
    (*B)(3,0) = h*gain_cmdvelocity*kp_torque;
  }

  {
    //Do the simulation:
    int nofsteps = round(h/(this->h));//Nosteps 
    if(m_dynamicsWorld)
    {
      int numSimSteps = m_dynamicsWorld->stepSimulation(h,nofsteps,this->h);//No interpolation used
      btTransform &chassistr = m_carChassis->getWorldTransform();
      btVector3 &chassisorig = chassistr.getOrigin();
      btMatrix3x3 &chbasis = chassistr.getBasis();
      btVector3 rpybasis;
      chbasis.getEulerZYX(rpybasis[2],rpybasis[1],rpybasis[0]);
      //Set the state
      xb<<chassisorig.x(), chassisorig.z(), rpybasis[2], m_vehicle->getCurrentSpeedKmHour()*(18.0/5.0)*(1/6.4817);//The gain factor 6.4817 is from MATLAB

      /*if(numSimSteps == 0)
      {
        cout<<"No simulation done"<<endl;
        cout<<"h: "<<h<<"\t"<<nofsteps<<"\t"<<(this->h)<<endl;
      }
      */
      this->x = xb;
      this->t += nofsteps*(this->h);
      //std::cout<<"this->x: "<<(this->x).transpose()<<"\t Time: "<<(this->t)<<endl;
      //cout<<(this->x).transpose()<<"\t"<<u.transpose()<<"\t"<<(this->t)<<endl;
    }
  }
  
  return 1;
}

double Bulletrccar::Step3(Vector4d &xb, const Vector2d &u,
    const Vector4d &w, double h,
    const VectorXd *p,Matrix4d *A, Matrix42d *B, Matrix4pd *C, Matrix4d *D)
{
  //Convert w into forces that can be applied:
  m_carChassis->applyCentralForce(btVector3(w[0], w[1], w[2]));
  m_carChassis->applyTorque(btVector3(0, w[3], 0));//Can be all dims too or change this acc to needs
  double result = this->Step1(xb, u, h, p, A, B, C);
  return result;
}

bool Bulletrccar::reset(const Vector4d &x, double t)
{
  if(m_dynamicsWorld)
  {
    gVehicleSteering = 0.f;
    gEngineForce = 0.f;
    gBreakingForce = 0.f;

    //Set to the specified state:
    btTransform cartransform(btQuaternion(x[3],0,0),btVector3(x[0],0,x[2]));
    m_carChassis->setCenterOfMassTransform(cartransform);
    btVector3 carvel(x[3]*6.4817,0,0);
    carvel = cartransform*carvel;
    m_carChassis->setLinearVelocity(carvel);
    m_carChassis->setAngularVelocity(btVector3(0,0,0));

    m_dynamicsWorld->getBroadphase()->getOverlappingPairCache()->cleanProxyFromPairs(m_carChassis->getBroadphaseHandle(),m_dynamicsWorld->getDispatcher());
    m_vehicle->resetSuspension();
    //Can synchronize wheels with the new transform or ignore
  }
  this->x = x;
  this->t = t;
}

bool Bulletrccar::NoiseMatrix(Matrix4d &Q, double t, const Vector4d &x, const Vector2d &u,
    double dt, const VectorXd *p) {
  Q.setZero();
  return true;
}

