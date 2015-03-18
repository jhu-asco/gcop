#include "bulletrccar.h"
#include <limits>

using namespace gcop;
using namespace Eigen;

Bulletrccar::Bulletrccar(BulletWorld&		m_world_, vector<double> *zs_) : 
  Rccar(3),m_carChassis(0), chassisShape(0), m_vehicle(0), m_wheelShape(0)
  ,gEngineForce(0),gBreakingForce(0.0),maxEngineForce(100.0),maxBreakingForce(10.0)
  ,gVehicleSteering(0.0),steeringClamp(0.5), velocityClamp(5.0)
  //,carmass(20.0), wheelRadius(0.07), wheelWidth(0.06f),wheelFriction(std::numeric_limits<double>::max())
  ,carmass(20.0), wheelRadius(0.07), wheelWidth(0.06f),wheelFriction(1.0)
  ,suspensionStiffness(250.0),suspensionDamping(100.0),suspensionCompression(0.02)
  ,rollInfluence(0.1),suspensionRestLength(0.122)
  ,m_defaultContactProcessingThreshold(1e5)
  , gain_cmdvelocity(1), kp_torque(25), kp_steer(0.1), initialz(0.15)
  ,m_world(m_world_), zs(zs_), reset_drivevel(false), initialstate(0)
{
//,rightIndex(0),upIndex(1),forwardIndex(2)
  if(!m_world.IsZupAxis())
  {
    //If world up axis is y then set the cars up axis also y
    rightIndex = 0;
    upIndex = 1;
    forwardIndex = 2;
    wheelDirectionCS0.setValue(0,-1,0);
    wheelAxleCS.setValue(-1,0,0);
    this->car_halfdims[0] = 0.19;//width
    this->car_halfdims[1] = 0.1;//height
    this->car_halfdims[2] =  0.32;//length
    this->l = 2*(this->car_halfdims[2]-wheelRadius); 

    //Set offset transform between car in GCOP coordsys and bullet
    offsettrans.setOrigin(btVector3(0,0,0));
    offsettrans.setRotation(btQuaternion(0,0,0));
  }
  else
  {
    //If the world up axis is z set the car z axis up forward axis as x and right axis as y
    forwardIndex = 1;
    rightIndex = 0;
    upIndex = 2;
    wheelDirectionCS0.setValue(0,0,-1);
    wheelAxleCS.setValue(1,0,0);//Can change on testing
    this->car_halfdims[0] = 0.19;//width
    this->car_halfdims[1] = 0.32;//length
    this->car_halfdims[2] =  0.1;//height
    this->l = 2*(this->car_halfdims[0]-wheelRadius); 

    //Set offset transform between car in GCOP coordsys and bullet
    offsettrans.setOrigin(btVector3(0,0,0));
    offsettrans.setRotation(btQuaternion(0,0,M_PI/2));
    //offsettrans.setRotation(btQuaternion(0,0,0));
    //Print Offset Trans:
    btMatrix3x3 basis = offsettrans.getBasis();
    offsettransinv = offsettrans.inverse();
  }
  //Load Car chassis and ray caster vehicle:
  {
    chassisShape = new btBoxShape(btVector3(car_halfdims[0],car_halfdims[1],car_halfdims[2]));//Need to use half extents of the box to create the chassis shape  in BULLET !!	
    //m_world.m_collisionShapes.push_back(chassisShape);

    btTransform initialtr;//temporary transform
    initialtr.setIdentity();
    assert(chassisShape);
    m_carChassis = m_world.LocalCreateRigidBody(carmass,initialtr,chassisShape);//(chassis Rigid Body)

    //Loading Wheel Shapes:
    m_wheelShape = new btCylinderShapeX(btVector3(wheelWidth,wheelRadius,wheelRadius));

    // Ray casting helper class for the vehicle
    m_vehicleRayCaster = new btDefaultVehicleRaycaster(m_world.m_dynamicsWorld);
    // Ray Casting Vehicle is made from existing rigid body. It converts a rigid body into a vehicle by attaching four rays looking downwards at four ends which are equivalent to wheels and uses the forces generated to move the rigid body.
    m_vehicle = new btRaycastVehicle(m_tuning,m_carChassis,m_vehicleRayCaster);

    ///never deactivate the vehicle from computing collisions etc. Usually the bodies are deactivated after some time if they are not moving
    m_carChassis->setActivationState(DISABLE_DEACTIVATION);

    (m_world.m_dynamicsWorld)->addVehicle(m_vehicle);
  }

  /// create wheel connections
  if(!m_world.IsZupAxis())
  {
    double connectionHeight = 0.17 - 0.03 - car_halfdims[1];// Suspension tip connection height - ground clearance - half the height of car


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
  else
  {
    double connectionHeight = 0.17 - 0.03 - car_halfdims[2];// Suspension tip connection height - ground clearance - half the height of car


    bool isFrontWheel=true;
    //choose coordinate system
    m_vehicle->setCoordinateSystem(rightIndex,upIndex,forwardIndex);

    //Front Right Wheel
    btVector3 connectionPointCS0(0.19,car_halfdims[1]-wheelRadius,connectionHeight);

    m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

    //Front Left Wheel
    connectionPointCS0 = btVector3(-0.19,car_halfdims[1]-wheelRadius,connectionHeight);
    m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

    //Back Right Wheel
    connectionPointCS0 = btVector3(0.19,-car_halfdims[1]+wheelRadius,connectionHeight);

    isFrontWheel = false;
    m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

    //Back Left Wheel
    connectionPointCS0 = btVector3(-0.19,-car_halfdims[1]+wheelRadius,connectionHeight);
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

  U.bnd = true;
  U.lb[0] = -velocityClamp;
  U.ub[0] = velocityClamp;//Max cmd
  U.lb[1] = -steeringClamp;
  U.ub[1] = steeringClamp;
}

double Bulletrccar::Step(Vector4d &xb, double t, const Vector4d &xa,
    const Vector2d &u, double h, const VectorXd *p,
    Matrix4d *A, Matrix42d *B, Matrix4pd *C)
{
//  cout<<"Step Called: "<<endl;
  this->reset(xa,t);
  return this->Step_internalinput(xb,u,h,p,A,B,C);
}

double Bulletrccar::Step_internalinput(Vector4d& xb, const Vector2d& u, 
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
 
  /*Set Linearization parameters:
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
    (*A)(3,3) = -h*kp_torque;
  }


  if (B) {
    B->setZero();
    (*B)(2,1) =  h*kp_steer*v/(l*pow(cos(gVehicleSteering),2));
    (*B)(3,0) = h*gain_cmdvelocity*kp_torque;
  }
  */

  {
    //Do the simulation:
    int nofsteps = round(h/(this->h));//Nosteps 
    for(int count_steps = 0;count_steps < nofsteps; count_steps++)
    {
      //Set Controls
      //double vel = m_carChassis->getLinearVelocity().length();
      double vel = m_vehicle->getCurrentSpeedKmHour()*(1.0/3.6);
      if(reset_drivevel)
      {
        vel = m_carChassis->getLinearVelocity().length();
        //cout<<"Vel:(Chassis vs Vehicle vel) "<<vel<<"\t"<<m_vehicle->getCurrentSpeedKmHour()*(1.0/3.6)<<endl;//#DEBUG
        reset_drivevel = false;
      }
      gEngineForce = kp_torque*(gain_cmdvelocity*u[0] - vel);
      //cout<<"More DEBUG: "<<gain_cmdvelocity<<"\t"<<u[0]<<"this_vel: "<<vel<<"\t"<<kp_torque<<"\t"<<gEngineForce<<endl;//#DEBUG
      gVehicleSteering = (1-kp_steer)*gVehicleSteering + kp_steer*u[1];
      gVehicleSteering = gVehicleSteering>steeringClamp?steeringClamp:(gVehicleSteering<-steeringClamp)?(-steeringClamp):gVehicleSteering;//Clamp Vehicle Steering based on bound
      //std::cout<<"Engine Force: "<<gEngineForce<<"\tgVehicleSteering: "<<gVehicleSteering<<std::endl;
      if(!m_vehicle)
      {
        cout<<"M_Vehicle not defined"<<endl;
      }
      int wheelIndex = 2;
      m_vehicle->applyEngineForce(gEngineForce,wheelIndex);

      wheelIndex = 3;
      m_vehicle->applyEngineForce(gEngineForce,wheelIndex);
      wheelIndex = 0;
      m_vehicle->setSteeringValue(gVehicleSteering,wheelIndex);
      wheelIndex = 1;
      m_vehicle->setSteeringValue(gVehicleSteering,wheelIndex);
      (m_world.m_dynamicsWorld)->stepSimulation((this->h),1,this->h);//No interpolation used Doing control of car at 50 Hz Depends on the physics engine timestep 100Hz
    }

    btTransform chassistr = m_carChassis->getWorldTransform();
    if(m_world.IsZupAxis())
    {
      chassistr = offsettransinv*chassistr*offsettrans;
    }
    btVector3 &chassisorig = chassistr.getOrigin();
    btMatrix3x3 &chbasis = chassistr.getBasis();
    btVector3 rpybasis;
    chbasis.getEulerZYX(rpybasis[2],rpybasis[1],rpybasis[0]);
    //Set the state
    if(!m_world.IsZupAxis())
    {
      xb<<chassisorig.x(), chassisorig.z(), rpybasis[1], m_vehicle->getCurrentSpeedKmHour()*(1.0/3.6);//The gain factor 1 is from MATLAB (6.4817)
    }
    else
    {
      xb<<chassisorig.x(), chassisorig.y(), rpybasis[2], m_vehicle->getCurrentSpeedKmHour()*(1.0/3.6);//The gain factor 1 is from MATLAB
    } 
    //#DEBUG: USEFUL
    //cout<<"Before Step:"<<endl;
    //cout<<(this->x).transpose()<<"\t"<<u.transpose()<<"\t"<<(this->t)<<endl;
    //btVector3 chassisvel = m_carChassis->getLinearVelocity();
    //cout<<"Compare vehicle Vel with chassis vel: "<<m_vehicle->getCurrentSpeedKmHour()*(18.0/5.0)*(1.0/12.9634)<<"\t"<<chassisvel.length()<<endl;
    //cout<<"Velocity: "<<chassisvel.x()<<"\t"<<chassisvel.y()<<"\t"<<chassisvel.z()<<endl;
    /*if(numSimSteps == 0)
      {
      cout<<"No simulation done"<<endl;
      cout<<"h: "<<h<<"\t"<<nofsteps<<"\t"<<(this->h)<<endl;
      }
     */
    this->x = xb;
    this->t += nofsteps*(this->h);
    //Set internal zs:
    if(this->zs)
    {
      count_zs++;
      assert(count_zs < zs->size());//#ASSERT
      (*zs)[count_zs] = chassisorig.z();
    }
    //cout<<"After Step:"<<endl;
    //std::cout<<(this->x).transpose()<<endl;//#DEBUG
    //std::cout<<chassisorig.x()<<"\t"<<chassisorig.y()<<"\t"<<chassisorig.z()<<std::endl;
  }
  
  return 1;
}

double Bulletrccar::Step_noise(Vector4d &xb, const Vector2d &u,
    const Vector4d &w, double h,
    const VectorXd *p,Matrix4d *A, Matrix42d *B, Matrix4pd *C, Matrix4d *D)
{
  //Convert w into forces that can be applied[NOT WORKING]:
  m_carChassis->applyCentralForce(btVector3(w[0], w[1], w[2]));
  m_carChassis->applyTorque(btVector3(0, w[3], 0));//Can be all dims too or change this acc to needs
  double result = this->Step_internalinput(xb, u, h, p, A, B, C);
  return result;
}

//Finds the closest car pose corresponding to the given pose
void Bulletrccar::setinitialstate(const Bulletrccar::CarState &inputstate, Vector4d &x)
{
  if(!initialstate)
    initialstate = new Bulletrccar::CarState;
  //Copy initialstate to inputstate:
  (*initialstate) = inputstate;
  {
    btTransform temp;
    temp.setIdentity();
    m_carChassis->setCenterOfMassTransform(temp);
  }
  btMatrix3x3 chbasis = inputstate.cartransform.getBasis();
  const btVector3 &source = inputstate.cartransform.getOrigin();// + (-0.2)*chbasis[2];
  btVector3 &carorigin = initialstate->cartransform.getOrigin();
  //const btVector3 &target = source + (-1.0)*chbasis[2];//-ve z dirxn
  const btVector3 &target = source + (-1.0)*btVector3(0,0,2.0);//-ve z dirxn
  btVector3 rpybasis;
  chbasis.getEulerZYX(rpybasis[2],rpybasis[1],rpybasis[0]);
  //#DEBUG:
  {
    cout<<"source: "<<source.x()<<"\t"<<source.y()<<"\t"<<source.z()<<"\t"<<endl;
    cout<<"target: "<<target.x()<<"\t"<<target.y()<<"\t"<<target.z()<<"\t"<<endl;
  }
  //do a ray test:
  btCollisionWorld::ClosestRayResultCallback rayCallback(source,target);
  (m_world.m_dynamicsWorld)->rayTest(source, target, rayCallback);
  if(rayCallback.hasHit())
  {
    cout<<"Distance to terrain: "<<(rayCallback.m_closestHitFraction)<<endl;
    //carorigin = carorigin + (suspensionRestLength + wheelRadius - 0.05 - 0.5*rayCallback.m_closestHitFraction)*chbasis[2];//Move the  car to intersection betn car and terrain + suspensionRestLength+wheelRadius - some nominal height
    carorigin = rayCallback.m_hitPointWorld + (suspensionRestLength + wheelRadius)*chbasis[2];//Move the  car to intersection betn car and terrain + suspensionRestLength+wheelRadius - some nominal height
    cout<<"hitPoint: "<<(rayCallback.m_hitPointWorld).x()<<"\t"<<(rayCallback.m_hitPointWorld).y()<<"\t"<<(rayCallback.m_hitPointWorld).z()<<"\t"<<endl;
    cout<<"chbasis: "<<(chbasis[2]).x()<<"\t"<<(chbasis[2]).y()<<"\t"<<(chbasis[2]).z()<<"\t"<<endl;
    cout<<"new origin: "<<carorigin[0]<<"\t"<<carorigin[1]<<"\t"<<carorigin[2]<<"\t"<<endl;
  }
  //verify all wheels are in contact:
  m_carChassis->setCenterOfMassTransform(initialstate->cartransform);
  bool allwheels_incontact = false;
  for(int count_wheels = 0;count_wheels < 4; count_wheels++)
  {
    btWheelInfo &wheelinfo = m_vehicle->m_wheelInfo[count_wheels];
    m_vehicle->rayCast(wheelinfo);
    cout<<count_wheels<<" wheel contact test: "<<wheelinfo.m_raycastInfo.m_isInContact<<endl;
    if(!wheelinfo.m_raycastInfo.m_isInContact)
    {
      allwheels_incontact = false;
      break;
    }
    //m_vehicle->updateWheelTransform(count_wheels, false);
  }

  //set the state:
  if(!m_world.IsZupAxis())
  {
    x<<carorigin.x(), carorigin.z(), rpybasis[1], initialstate->carlinearvel.length();
  }
  else
  {
    x<<carorigin.x(), carorigin.y(), rpybasis[2], initialstate->carlinearvel.length();
  }
  //set internal zs[0]:
  if(zs)
    (*zs)[0] = carorigin.z();
  
  //#DEBUG:
  {
  //  cout<<"carlinearvel: "<<(initialstate->carlinearvel).x()<<"\t"<<(initialstate->carlinearvel).y()<<"\t"<<(initialstate->carlinearvel).z()<<"\t"<<endl;
    //cout<<"carangularvel: "<<(initialstate->carangularvel).x()<<"\t"<<(initialstate->carangularvel).y()<<"\t"<<(initialstate->carangularvel).z()<<"\t"<<endl;
    cout<<"State: "<<x.transpose()<<endl;
  }
}
void Bulletrccar::setinitialstate(Vector4d &x)
{
  if(!initialstate)
    initialstate = new Bulletrccar::CarState;
  initialstate->cartransform = m_carChassis->getWorldTransform();
  initialstate->carlinearvel = m_carChassis->getLinearVelocity();
  initialstate->carangularvel = m_carChassis->getAngularVelocity();

  btTransform chassistr = initialstate->cartransform;
  if(m_world.IsZupAxis())
  {
    chassistr = offsettransinv*chassistr*offsettrans;
  }
  btVector3 &chassisorig = chassistr.getOrigin();
  btMatrix3x3 &chbasis = chassistr.getBasis();
  btVector3 rpybasis;
  chbasis.getEulerZYX(rpybasis[2],rpybasis[1],rpybasis[0]);
  //Set the state
  if(!m_world.IsZupAxis())
  {
    x<<chassisorig.x(), chassisorig.z(), rpybasis[1], initialstate->carlinearvel.length();
  }
  else
  {
    x<<chassisorig.x(), chassisorig.y(), rpybasis[2], initialstate->carlinearvel.length();
  }
  //set internal zs[0]:
  if(zs)
    (*zs)[0] = chassisorig.z();
  
  //#DEBUG:
  {
    cout<<"carlinearvel: "<<(initialstate->carlinearvel).x()<<"\t"<<(initialstate->carlinearvel).y()<<"\t"<<(initialstate->carlinearvel).z()<<"\t"<<endl;
    cout<<"carangularvel: "<<(initialstate->carangularvel).x()<<"\t"<<(initialstate->carangularvel).y()<<"\t"<<(initialstate->carangularvel).z()<<"\t"<<endl;
    cout<<"State: "<<x.transpose()<<endl;
  }
}

bool Bulletrccar::reset(const Vector4d &x, double t)
{
  if(m_world.m_dynamicsWorld)
  {
    gVehicleSteering = 0;
    gEngineForce = 0;
    gBreakingForce = 0;

    // Delete the Vehicle and re initialize:
    {
      m_world.m_dynamicsWorld->removeVehicle(m_vehicle);
      delete m_vehicle;
      m_world.m_dynamicsWorld->removeRigidBody(m_carChassis);
      delete m_carChassis;


      btTransform initialtr;//temporary transform
      initialtr.setIdentity();

      m_carChassis = m_world.LocalCreateRigidBody(carmass,initialtr,chassisShape);//(chassis Rigid Body)

      // Ray casting helper class for the vehicle
      //m_vehicleRayCaster = new btDefaultVehicleRaycaster(m_world.m_dynamicsWorld);

      m_vehicle = new btRaycastVehicle(m_tuning,m_carChassis,m_vehicleRayCaster);

      ///never deactivate the vehicle from computing collisions etc. Usually the bodies are deactivated after some time if they are not moving
      m_carChassis->setActivationState(DISABLE_DEACTIVATION);

      (m_world.m_dynamicsWorld)->addVehicle(m_vehicle);

      /// create wheel connections
      if(!m_world.IsZupAxis())
      {
        double connectionHeight = 0.17 - 0.03 - car_halfdims[1];// Suspension tip connection height - ground clearance - half the height of car


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
      else
      {
        double connectionHeight = 0.17 - 0.03 - car_halfdims[2];// Suspension tip connection height - ground clearance - half the height of car


        bool isFrontWheel=true;
        //choose coordinate system
        m_vehicle->setCoordinateSystem(rightIndex,upIndex,forwardIndex);

        //Front Right Wheel
        btVector3 connectionPointCS0(0.19,car_halfdims[1]-wheelRadius,connectionHeight);

        m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

        //Front Left Wheel
        connectionPointCS0 = btVector3(-0.19,car_halfdims[1]-wheelRadius,connectionHeight);
        m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

        //Back Right Wheel
        connectionPointCS0 = btVector3(0.19,-car_halfdims[1]+wheelRadius,connectionHeight);

        isFrontWheel = false;
        m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

        //Back Left Wheel
        connectionPointCS0 = btVector3(-0.19,-car_halfdims[1]+wheelRadius,connectionHeight);
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

    if(initialstate)
    {
      m_carChassis->setCenterOfMassTransform(initialstate->cartransform);
      m_carChassis->setLinearVelocity(initialstate->carlinearvel);
      m_carChassis->setAngularVelocity(initialstate->carangularvel);
    }
    else
    {
      //Set to the specified state:
      if(!m_world.IsZupAxis())
      {
        btTransform cartransform(btQuaternion(x[2],0,0),btVector3(x[0], initialz,x[1]));
        m_carChassis->setCenterOfMassTransform(cartransform);
        double v = x[3];
        //double v = x[3];
        btVector3 carvel(-v*sin(x[2]),0,v*cos(x[2]));
        m_carChassis->setLinearVelocity(carvel);
      }
      else
      {
        btTransform cartransform(btQuaternion(0,0,x[2]),btVector3(x[0], x[1], initialz));
        //cartransform = offsettrans*cartransform*offsettransinv;
        m_carChassis->setCenterOfMassTransform(cartransform);
        double v = x[3];
        //double v = x[3];
        //btVector3 carvel(-v*sin(x[2]),v*cos(x[2]),0);
        btVector3 carvel(0,0,0);
        m_carChassis->setLinearVelocity(carvel);
        //      cout<<"Car Vel: "<<carvel.x()<<"\t"<<carvel.y()<<"\t"<<carvel.z()<<endl;
      }
      m_carChassis->setAngularVelocity(btVector3(0,0,0));
    }

    m_world.Reset();
    //(m_world.m_dynamicsWorld)->getBroadphase()->getOverlappingPairCache()->cleanProxyFromPairs(m_carChassis->getBroadphaseHandle(),(m_world.m_dynamicsWorld)->getDispatcher());
    //if(m_vehicle)
    //{
    // m_vehicle->resetSuspension();
    //Can synchronize wheels with the new transform
    /*for (int i=0;i<m_vehicle->getNumWheels();i++)
      {
    //synchronize the wheels with the (interpolated) chassis worldtransform
    m_vehicle->updateWheelTransform(i,true);
    }
     */
    //m_vehicle->updateVehicle(0);
    //}
    //set reset_drivevel to true:
    reset_drivevel = true;
  }
  this->x = x;
  this->t = t;
  //Set internal zs also:
  if(zs)
  {
    if(!initialstate)
      (*zs)[0] = initialz;
    count_zs = 0;
  }
}

bool Bulletrccar::NoiseMatrix(Matrix4d &Q, double t, const Vector4d &x, const Vector2d &u,
    double dt, const VectorXd *p) {
  Q.setZero();
  return true;
}

