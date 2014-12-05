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
#include <LinearMath/btDefaultMotionState.h>
#include <BulletDynamics/Dynamics/btDynamicsWorld.h>
#include <BulletCollision/CollisionShapes/btCollisionShape.h>
#include <BulletCollision/CollisionShapes/btCollisionShape.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <BulletDynamics/ConstraintSolver/btConstraintSolver.h>

#include <iostream>
#include <stdio.h>

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

      //Creates a triangular mesh from vertex data. The number of indices should be equal to 3 times the number of triangles
      // The vertex data can be smaller than 3*noftriangles since different triangles can share vertices
      // By default assumes each index corresponds to a vertex
      btCollisionShape *CreateMeshFromData(float *verts, int *inds, int noftriangles, int nofverts = (3*noftriangles))
      {
        if(!inds || !verts)
          return;
        int vertStride = sizeof(float); 
        int indexStride = 3*sizeof(int);
        m_indexVertexArrays = new btTriangleIndexVertexArray(noftriangles,
            inds,
            indexStride,
            nofverts,groundverts,vertStride);
        cout<<"2"<<endl;

        bool useQuantizedAabbCompression = true;
        btCollisionShape *groundShape = new btBvhTriangleMeshShape(m_TriMesh,useQuantizedAabbCompression);
        return groundShape;
      }

      //Loads a mesh from Binary STL Files. Modified from LoadMeshFromSTL used by Bullet3 Demo
			btCollisionShape *CreateMeshFromSTL(const char *filename)
      {
        FILE* file = fopen(filename,"rb");
        if(file)
        {
          int size=0;
          if (fseek(file, 0, SEEK_END) || (size = ftell(file)) == EOF || fseek(file, 0, SEEK_SET))
          {
            printf("Error: Cannot access file to determine size of %s\n", filename);
          } 
          else
          {
            if (size)
            {
              printf("Open STL file of %d bytes\n",size);
              char* memoryBuffer = new char[size+1];
              int actualBytesRead = fread(memoryBuffer,1,size,file);
              if (actualBytesRead!=size)
              {
                printf("Error reading from file %s",filename);
              } 
              else
              {
                int numTriangles = *(int*)&memoryBuffer[80];

                if (numTriangles)
                {
                  {
                    //perform a sanity check instead of crashing on invalid triangles/STL files
                    int expectedBinaryFileSize = numTriangles* 50 + 84;
                    if (expectedBinaryFileSize != size)
                    {
                      return 0;
                    }

                  }
                  float *verts = (float *)malloc(3*numTriangles*sizeof(float));//No compression used Can use trees if needed later
                  int *inds = (int *)malloc(3*numTriangles*sizeof(int));
                  int i;
#pragma omp parallel for private(i)
                  for (i=0;i<numTriangles;i++)
                  {
#pragma omp section
                    char* curPtr = &memoryBuffer[84+i*50];
                    memcpy(verts+3*i,curPtr+12,36);
#pragma omp section
                    inds[3*i+0] = i;
                    inds[3*i+1] = i+1;
                    inds[3*i+2] = i+2;
                  }
                }
                return CreateMeshFromData(verts, inds, numTriangles);

                delete[] memoryBuffer;
              }
            }
            fclose(file);
          }
        }
        return 0;
      }

			//Simulates the time dt in seconds with the specified number of substeps
			void Step(float dt, int substeps)
			{
				if(!m_dynamicsWorld)
        {
          cerr<<"m_dynamicsWorld not defined"<<endl;
          return;
        }
				float fixedstepsize = dt/substeps;
				int numSimSteps = m_dynamicsWorld->stepSimulation(dt,substeps,fixedstepsize);//No interpolation used
			}

			//Reset Bullet Physics Engine and clear collision info and reset constraint solver
			void Reset()
      {
        if(!m_dynamicsWorld)
				{
					cerr<<"m_dynamicsWorld not defined"<<endl;
					return;
				}
        m_overlappingPairCache->resetPool(m_dispatcher);
				m_constraintSolver->reset();
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

				delete m_dynamicsWorld;

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
