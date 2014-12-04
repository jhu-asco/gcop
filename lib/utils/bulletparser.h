/** Bullet Helper Class:
 * This class is derived from the Bullet Demo application to provide 
 * helper functions to setup a default world (The code here is just copy paste from above sources without addition OpenGL stuff).
 * Once a world is setup, you can load different models into it and perform iterations.
 *
 * Author: Gowtham Garimella (ggarime1@jhu.edu)
 * 
 */
#include <BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <btBulletDynamicsCommon.h>
#include <BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h>
#include <LinearMath/btQuickprof.h>
#include <LinearMath/btDefaultMotionState.h>
#include <LinearMath/btSerializer.h>
#include <LinearMath/btIDebugDraw.h>
#include <BulletDynamics/Dynamics/btDynamicsWorld.h>

#include <BulletDynamics/ConstraintSolver/btPoint2PointConstraint.h>//picking
#include <BulletDynamics/ConstraintSolver/btGeneric6DofConstraint.h>//picking

#include <BulletCollision/CollisionShapes/btCollisionShape.h>
#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/CollisionShapes/btSphereShape.h>
#include <BulletCollision/CollisionShapes/btCompoundShape.h>
#include <BulletCollision/CollisionShapes/btUniformScalingShape.h>
#include <BulletDynamics/ConstraintSolver/btConstraintSolver.h>

#include <iostream>



using namespace std;

namespace BulletWorld{
	class BulletWorld
	{
		protected:
			btAlignedObjectArray<btCollisionShape*> m_collisionShapes;///<Stores all the collision shapes in a world for drawing and other purposes

			class btBroadphaseInterface*	m_overlappingPairCache;///<Broadphase collision checker

			class btCollisionDispatcher*	m_dispatcher;///<Collision checker

			class btConstraintSolver*	m_constraintSolver;///<Constraint Solver

			class btDefaultCollisionConfiguration* m_collisionConfiguration;///<Default Collision configuration used from bullet demos

			btDynamicsWorld*		m_dynamicsWorld;///<this is the most important class It represents the Discrete Dynamics World 

      btVector3 worldMin;
      btVector3 worldMax;

		public:
      BulletWorld(btVector3 worldmin_ = btVector3(-1000,-1000,-1000), btVector3 worldMax_ = btVector3(1000,1000,1000), btVector3 gravity_ = btVector3(0,0,-9.81)):
        m_collisionConfiguration(new btDefaultCollisionConfiguration())
        ,m_dispatcher(new btCollisionDispatcher(m_collisionConfiguration))
        ,worldMin(worldMin_), worldMax(worldMax_)
        ,m_overlappingPairCache(new btAxisSweep3(worldMin,worldMax))
        ,m_constraintSolver(new btSequentialImpulseConstraintSolver())
        ,m_dynamicsWorld(new btDiscreteDynamicsWorld(m_dispatcher,m_overlappingPairCache,m_constraintSolver,m_collisionConfiguration))
			{
				m_dynamicsWorld->setGravity(gravity_);//Sets Gravity Vector
			}

			btRigidBody*  localCreateRigidBody(float mass, const btTransform& startTransform,btCollisionShape* shape, bool use_motion_state = false)
			{
				btAssert((!shape || shape->getShapeType() != INVALID_SHAPE_PROXYTYPE));

				//rigidbody is dynamic if and only if mass is non zero, otherwise static
				bool isDynamic = (mass != 0.f);


				btVector3 localInertia(0,0,0);
				if (isDynamic)
					shape->calculateLocalInertia(mass,localInertia);

				//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
        if(use_motion_state)
        {
          btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);

          btRigidBody::btRigidBodyConstructionInfo cInfo(mass,myMotionState,shape,localInertia);

          btRigidBody* body = new btRigidBody(cInfo);
          body->setContactProcessingThreshold(m_defaultContactProcessingThreshold);
        }
        else
        {
          btRigidBody* body = new btRigidBody(mass,0,shape,localInertia);
          body->setWorldTransform(startTransform);
        }

				m_dynamicsWorld->addRigidBody(body);

				return body;
			}


      //Copied Code from BulletMeshShape.cc:
			void Mesh(float *groundverts, int *groundinds)
			{
        if(
 				m_TriMesh = new btTriangleMesh();//Create a Triangle Mesh
				ROS_INFO("2");




				//assert(submeshcount == 1);
				//float *groundvert1 = groundverts[0];
				//float *groundinds1 = groundinds[0];//1st submesh verts and inds
				int vertStride = 3*sizeof(float); 
				int indexStride = 3*sizeof(int);
				int numIndices = groundmesh->GetIndexCount();//This is not the total Triangles just for debugging
				int numVertices = groundmesh->GetVertexCount();
				ROS_INFO("Vertex count: %d",numVertices);
				ROS_INFO("Index Count: %d",numIndices);
				// Scale the vertex data
				for (unsigned int j = 0;  j < numVertices; j++)
				{
					groundverts[j*3+0] = groundverts[j*3+0] *
							_sdf->Get<gazebo::math::Vector3>("scale").x;
					groundverts[j*3+1] = groundverts[j*3+1] *
							_sdf->Get<gazebo::math::Vector3>("scale").y;
					groundverts[j*3+2] = groundverts[j*3+2] *
							_sdf->Get<gazebo::math::Vector3>("scale").z;
				}
				// Create the Bullet trimesh
				for (unsigned int j = 0; j < numIndices; j += 3)
				{
					btVector3 bv0(groundverts[groundinds[j]*3+0],
							groundverts[groundinds[j]*3+1],
							groundverts[groundinds[j]*3+2]);

					btVector3 bv1(groundverts[groundinds[j+1]*3+0],
							groundverts[groundinds[j+1]*3+1],
							groundverts[groundinds[j+1]*3+2]);

					btVector3 bv2(groundverts[groundinds[j+2]*3+0],
							groundverts[groundinds[j+2]*3+1],
							groundverts[groundinds[j+2]*3+2]);

					m_TriMesh->addTriangle(bv0, bv1, bv2);
					/*gzdbg<<bv0.getX()<<"\t"<<bv0.getY()<<"\t"<<bv0.getZ()<<endl;
					gzdbg<<bv1.getX()<<"\t"<<bv1.getY()<<"\t"<<bv1.getZ()<<endl;
					gzdbg<<bv2.getX()<<"\t"<<bv2.getY()<<"\t"<<bv2.getZ()<<endl;
					*/
				}
				bool useQuantizedAabbCompression = true;
				//btCollisionShape *groundShape = new btBvhTriangleMeshShape(m_TriMesh,useQuantizedAabbCompression);
			 /* m_indexVertexArrays = new btTriangleIndexVertexArray(totalTriangles,
						groundinds,
						indexStride,
						totalVerts,groundverts,vertStride);
			 */
				btCollisionShape *groundShape = new btBvhTriangleMeshShape(m_TriMesh,useQuantizedAabbCompression);
				//btConvexShape* convexShape = new btConvexTriangleMeshShape(m_TriMesh, true);//Convex Shape does not make sense from my understanding as the meshes are not Convex in general
				//convexShape->setMargin(0.001f);
				//bParent->SetCollisionShape(convexShape);
				m_collisionShapes.push_back(groundShape);
				btTransform tr;
				tr.setIdentity();
				//tr.setOrigin(btVector3(0,0,-1.0f));
				localCreateRigidBody(0,tr,groundShape);//mass 0 i.e  immovable Object
				delete [] groundverts;
				delete [] groundinds;
			}

/*
			void LoadGround(const boost::shared_ptr<gazebo::physics::Model> modelptr)//Not going to use the ground plane data from model
			{
				//For now assume ground is a long plane of 50 x 50 m and 1 m box just to define a box in bullet
				btCollisionShape* groundShape = new btBoxShape(btVector3(25,25,0.5));
				m_collisionShapes.push_back(groundShape);
				float groundmass = 0;//immovable object
				btTransform tr;
				tr.setIdentity();
				tr.setOrigin(btVector3(0,0,-0.5f));//So the plane will have exactly 0m height
				btQuaternion qtr;
				//qtr.setEulerZYX(M_PI/2,0,0);
				qtr.setEulerZYX(0,0,0);
				tr.setRotation(qtr);
				localCreateRigidBody(0,tr,groundShape);//mass 0 i.e  immovable Object
			}
			*/
			void LoadCar(const boost::shared_ptr<gazebo::physics::Model> modelptr, btTransform initialtr_)//Not going to use the car data for now
			{
				//The car details are hardcoded
				float carmass = 120.f;//120 Kg
				btCollisionShape* chassisShape1 = new btBoxShape(btVector3(0.2f,0.615f, 0.09f));//Need to use half extents of the box to create the chassis shape  in BULLET !!	
				//	btCollisionShape* chassisShape2 = new btBoxShape(btVector3(0.3f,0.72f, 0.18f));//The documentation says to avoid creating rigidbodies one inside another
				btCompoundShape* compound = new btCompoundShape();
				btTransform localTrans;
				localTrans.setIdentity();
				//localTrans effectively shifts the center of mass with respect to the chassis
				localTrans.setOrigin(btVector3(0,0,0.2932));
				compound->addChildShape(localTrans,chassisShape1);
				//	compound->addChildShape(localTrans,chassisShape2);

        btTransform initialtr = initialtr_*offsettrans.inverse();

		    m_carChassis = localCreateRigidBody(carmass,initialtr,compound);//(chassis Rigid Body)
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

					float connectionHeight = wheelRadius;//Where wheels are connected with respect to the chassis world origin


					bool isFrontWheel=true;
					//choose coordinate system
					m_vehicle->setCoordinateSystem(rightIndex,upIndex,forwardIndex);

					//Front Right Wheel
					btVector3 connectionPointCS0(0.3048,0.7366-0.3683, connectionHeight);

					m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

					//Front Left Wheel
					connectionPointCS0 = btVector3(-0.3048,0.7366-0.3683, connectionHeight);

					m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

					//Back Right Wheel
					connectionPointCS0 = btVector3(0.3048,-0.3683, connectionHeight);

					isFrontWheel = false;
					m_vehicle->addWheel(connectionPointCS0,wheelDirectionCS0,wheelAxleCS,suspensionRestLength,wheelRadius,m_tuning,isFrontWheel);

					//Back Left Wheel
					connectionPointCS0 = btVector3(-0.3048, -0.3683, connectionHeight);

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

			//Applies the Forces to the Vehicle
			void ApplyVehicleForces(float engineForce, float breakingForce, float Steeringforce)
			{
				if(!m_vehicle)
				{
					cout<<"M_Vehicle not defined"<<endl;
				}
				if(abs(engineForce) > 0 && abs(breakingForce) > 0)
				{
					cerr<<"Both Engine and Breaking forces are nonzero"<<endl;
					return;
				}
				int wheelIndex = 2;
				m_vehicle->applyEngineForce(engineForce,wheelIndex);
				m_vehicle->setBrake(breakingForce,wheelIndex);

				wheelIndex = 3;
				m_vehicle->applyEngineForce(engineForce,wheelIndex);
				m_vehicle->setBrake(breakingForce,wheelIndex);
				wheelIndex = 0;
				m_vehicle->setSteeringValue(Steeringforce,wheelIndex);
				wheelIndex = 1;
				m_vehicle->setSteeringValue(Steeringforce,wheelIndex);
			}

			//Simulates the time dt in seconds with the specified number of substeps
			//Provides the chassistransformation at the end as an output
			//Later will add the input transformation to start the simulation from  this
			//will allow us to find the gradient approximation at a point
			void Simulate(float dt, int substeps, std::vector<btTransform> &outtrans, float &bodyVel)
			{
				if(!m_dynamicsWorld)
					cerr<<"m_dynamicsWorld not defined"<<endl;
				float fixedstepsize = dt/substeps;
				int numSimSteps = m_dynamicsWorld->stepSimulation(dt,substeps,fixedstepsize);//No interpolation used
				assert(outtrans.capacity() >= 5);//It should be able to hold atleast 5 elements
				if(m_carChassis && m_vehicle)
				{
					outtrans[0] = m_carChassis->getWorldTransform()*offsettrans;
					for(int countwheel = 0;countwheel < 4;countwheel++)
						outtrans[countwheel+1] = m_vehicle->getWheelTransformWS(countwheel);
					//chassistrans = m_carChassis->getWorldTransform();
					//printf("world pos = %f,%f,%f\n",float(chassistrans.getOrigin().getX()),float(chassistrans.getOrigin().getY()),float(chassistrans.getOrigin().getZ()));
					bodyVel = m_vehicle->getCurrentSpeedKmHour()*(18.0/5.0);//Convert to m/s
				}
			}
			//Reset Bullet Physics Engine and place the vehicle back at some given transformation
			void Reset(btTransform &transform, btVector3 &velocity, btVector3 &angularvelocity)
      {
				//There are atleast two ways of doing this.
				//One is to remove the rigidbody and add it back with the new state
				//You can also just try setWorldTransform and set the velocities to 0
				//Alternativey you can also do some contact clearing and stuff that can be seen from
				//DemoApplication: ClientResetScene
				//We are just trying to set the transform and velocities to zero at new position
				if(!m_dynamicsWorld || !m_carChassis)
				{
					cerr<<"the pointers dynamicsworld and carChassis not defined"<<endl;
					return;
				}
        btTransform initialtr = transform*offsettrans.inverse();
				m_carChassis->setCenterOfMassTransform(initialtr);//Reset body to initial Tr
				m_carChassis->setLinearVelocity(velocity);// Can take input to reset to a different pose
				m_carChassis->setAngularVelocity(angularvelocity);
        //if(m_dynamicsWorld->getBroadphase()->getOverlappingPairCache())
          //m_dynamicsWorld->getBroadphase()->getOverlappingPairCache()->cleanProxyFromPairs(m_carChassis->getBroadphaseHandle(),m_dispatcher);
				m_overlappingPairCache->resetPool(m_dispatcher);
				m_constraintSolver->reset();
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

      void Reset(btTransform &transform)
      {
        btVector3 zerovector(0,0,0);
			  this->Reset(transform, zerovector, zerovector);
      }
      

			~BulletWorld()//Destructor for Bullet Physics Engine
			{
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

				delete m_TriMesh;

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
	};
};
