#ifndef BULLET_WORLD
#define BULLET_WORLD
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
#include "BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h"
#include <BulletDynamics/ConstraintSolver/btConstraintSolver.h>

#include <iostream>
#include <stdio.h>

using namespace std;

namespace gcop{
	class BulletWorld
	{
		protected:
			class btBroadphaseInterface*	m_overlappingPairCache;///<Broadphase collision checker

			class btCollisionDispatcher*	m_dispatcher;///<Collision checker

			class btConstraintSolver*	m_constraintSolver;///<Constraint Solver

			class btDefaultCollisionConfiguration* m_collisionConfiguration;///<Default Collision configuration used from bullet demos


      btVector3 worldMin;
      btVector3 worldMax;
      bool usezupaxis;

    public: 
			btDynamicsWorld*		m_dynamicsWorld;///<this is the most important class It represents the Discrete Dynamics World 
      btScalar m_defaultContactProcessingThreshold;
			btAlignedObjectArray<btCollisionShape*> m_collisionShapes;///<Stores all the collision shapes in a world for drawing and other purposes


		public:
      BulletWorld(bool usezupaxis_ = false, btVector3 worldMin_ = btVector3(-1000,-1000,-1000), btVector3 worldMax_ = btVector3(1000,1000,1000)):
        m_defaultContactProcessingThreshold(BT_LARGE_FLOAT)
        ,worldMin(worldMin_), worldMax(worldMax_), usezupaxis(usezupaxis_)
      {
        m_collisionConfiguration = new btDefaultCollisionConfiguration();
        m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);
        m_overlappingPairCache = new btAxisSweep3(worldMin,worldMax);
        m_constraintSolver = new btSequentialImpulseConstraintSolver();
        m_dynamicsWorld = new btDiscreteDynamicsWorld(m_dispatcher,m_overlappingPairCache,m_constraintSolver,m_collisionConfiguration);
        if(usezupaxis)
        {
          m_dynamicsWorld->setGravity(btVector3(0,0,-9.81));//Sets Gravity Vector
          btVector3 gravity_ = m_dynamicsWorld->getGravity();
        }
        else
        {
          m_dynamicsWorld->setGravity(btVector3(0,-9.81,0));//Sets Gravity Vector
        }
      }

      bool IsZupAxis()
      {
        return usezupaxis;
      }

      void SetGravity(btVector3 gravity_)
      {
				m_dynamicsWorld->setGravity(gravity_);//Sets Gravity Vector
      }

			btRigidBody*  LocalCreateRigidBody(float mass, const btTransform& startTransform,btCollisionShape* shape, bool use_motion_state = false)
			{
				btAssert((!shape || shape->getShapeType() != INVALID_SHAPE_PROXYTYPE));

				//rigidbody is dynamic if and only if mass is non zero, otherwise static
				bool isDynamic = (mass != 0.f);


				btVector3 localInertia(0,0,0);
				if (isDynamic)
					shape->calculateLocalInertia(mass,localInertia);

        btRigidBody* body = NULL;
				//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
        if(use_motion_state)
        {
          btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);

          btRigidBody::btRigidBodyConstructionInfo cInfo(mass,myMotionState,shape,localInertia);

          body = new btRigidBody(cInfo);
          body->setContactProcessingThreshold(m_defaultContactProcessingThreshold);
        }
        else
        {
          body = new btRigidBody(mass,0,shape,localInertia);
          body->setWorldTransform(startTransform);
        }
        if(body != NULL)
          m_dynamicsWorld->addRigidBody(body);

				return body;
			}

      //Creates a triangular mesh from vertex data. The number of indices should be equal to 3 times the number of triangles
      // The vertex data can be smaller than 3*noftriangles since different triangles can share vertices
      // By default assumes each index corresponds to a vertex
      btCollisionShape *CreateMeshFromData(btScalar *verts, int *inds, int noftriangles, int nofverts)
      {
        if(!inds || !verts)
          return 0;
        int vertStride = 3*sizeof(btScalar); 
        int indexStride = 3*sizeof(int);
        btTriangleIndexVertexArray *m_indexVertexArrays = new btTriangleIndexVertexArray(noftriangles
                                                                                        ,inds
                                                                                        ,indexStride
                                                                                        ,nofverts,verts,vertStride);
        cout<<"2"<<endl;

        bool useQuantizedAabbCompression = true;
        btCollisionShape *collshape = new btBvhTriangleMeshShape(m_indexVertexArrays,useQuantizedAabbCompression);
        m_collisionShapes.push_back(collshape);
        return collshape;
      }

      //Loads a mesh from Binary STL Files. Modified from LoadMeshFromSTL used by Bullet3 Demo
			btCollisionShape *CreateMeshFromSTL(const char *filename, btVector3 scale = btVector3(1,1,1))
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
                  btScalar *verts = (btScalar *)malloc(3*3*numTriangles*sizeof(btScalar));//No compression used Can use trees if needed later 3 vertices per triangle, 3 indices per vertex
                  int *inds = (int *)malloc(3*numTriangles*sizeof(int));
                  int i;

                  {
                    float *vert_temp = (float *)malloc(3*3*numTriangles*sizeof(float));

#pragma omp parallel for private(i)
                    for (i=0;i<numTriangles;i++)
                    {
                      memcpy(&vert_temp[9*i],&memoryBuffer[96+i*50],36);//9 floats of 4 bytes each (3 loats per vertex)
                      inds[3*i+0] = 3*i;
                      inds[3*i+1] = 3*i+1;
                      inds[3*i+2] = 3*i+2;
                      //cout<<"Vertices["<<3*i<<"]: "<<verts[3*i]<<"\t"<<verts[3*i+1]<<"\t"<<verts[3*i+2]<<endl;
                      //cout<<"Vertices["<<(3*i+1)<<"]: "<<verts[3*(i+1)]<<"\t"<<verts[3*(i+1)+1]<<"\t"<<verts[3*(i+1)+2]<<endl;
                      //cout<<"Vertices["<<(3*i+2)<<"]: "<<verts[3*(i+2)]<<"\t"<<verts[3*(i+2)+1]<<"\t"<<verts[3*(i+2)+2]<<endl;
                    }

                    //Convert float vertices into  double vertices:
#pragma omp parallel for private(i)
                    for(i = 0; i < 9*numTriangles; i++)
                    {
                      verts[i] = double(vert_temp[i]);//explicit casting 
                    }
                  }

                  //Scale the vertices accordingly
                  if((scale - btVector3(1,1,1)).norm()>1e-6)
                  {
#pragma omp parallel for private(i)
                    for(i = 0;i < 9*numTriangles;i++)
                    {
                      verts[i] *= scale[i%3];
                    }
                  }

                  return CreateMeshFromData(verts, inds, numTriangles, 3*numTriangles);
                }
                delete[] memoryBuffer;
              }
            }
            fclose(file);
          }
        }
        return 0;
      }

			btCollisionShape *CreateHeightMap(btScalar length, btScalar width, const char *data, btScalar maxHeight = 100)
      {
        btHeightfieldTerrainShape* heightFieldShape;
        if(usezupaxis)
        {
          heightFieldShape = new btHeightfieldTerrainShape(width,length,data,maxHeight,2,false, false);
        }
        else
        {
          heightFieldShape = new btHeightfieldTerrainShape(width,length,data,maxHeight,1,false,false);
        }
        return (btCollisionShape*)heightFieldShape;
      }

			btCollisionShape *CreateGroundPlane(btScalar length, btScalar width, btScalar(*heightfunc)(btScalar, btScalar)=0,int subdivisions = 1)
      {
        cout<<"Subdivisions: "<<subdivisions<<endl;
        int i;

        const btScalar TRIANGLE_SIZE_X = (length/subdivisions);
        const btScalar TRIANGLE_SIZE_Y = (width/subdivisions);

        const int NUM_VERTS_X = subdivisions+1;
        const int NUM_VERTS_Y = subdivisions+1;
        const int totalVerts = NUM_VERTS_X*NUM_VERTS_Y;

        const int totalTriangles = 2*(NUM_VERTS_X-1)*(NUM_VERTS_Y-1);

        btScalar *m_vertices = new btScalar[3*totalVerts];
        int*	gIndices = new int[totalTriangles*3];

        cout<<"Number of Verts: "<<NUM_VERTS_X<<"\t"<<NUM_VERTS_Y<<endl;

#pragma omp parallel for private(i)
        for ( i=0;i<NUM_VERTS_X;i++)
        {
          for (int j=0;j<NUM_VERTS_Y;j++)
          {
            btScalar height = 0.f;//20.f*sinf(float(i)*wl)*cosf(float(j)*wl);
            if(heightfunc != 0)
            {
              height = heightfunc(i*TRIANGLE_SIZE_X, j*TRIANGLE_SIZE_Y);
            }

            m_vertices[3*(i+j*NUM_VERTS_X)] = (i-NUM_VERTS_X*0.5)*TRIANGLE_SIZE_X + 0.5*TRIANGLE_SIZE_X;
            if(usezupaxis)
            {
              m_vertices[3*(i+j*NUM_VERTS_X)+1] = (j-NUM_VERTS_Y*0.5)*TRIANGLE_SIZE_Y + 0.5*TRIANGLE_SIZE_Y;
              m_vertices[3*(i+j*NUM_VERTS_X)+2] = height;
            }
            else
            {
              m_vertices[3*(i+j*NUM_VERTS_X)+1] = height;
              m_vertices[3*(i+j*NUM_VERTS_X)+2] = (j-NUM_VERTS_Y*0.5)*TRIANGLE_SIZE_Y + 0.5*TRIANGLE_SIZE_Y;
            }
            //cout<<"Vertices: "<<3*(i+j*NUM_VERTS_X)<<"\t"<<m_vertices[3*(i+j*NUM_VERTS_X)]<<"\t"<<m_vertices[3*(i+j*NUM_VERTS_X)+1]<<"\t"<<m_vertices[3*(i+j*NUM_VERTS_X)+2]<<endl;
          }
        }

#pragma omp parallel for private(i)
        for ( i=0;i<NUM_VERTS_X-1;i++)
        {
          for (int j=0;j<NUM_VERTS_Y-1;j++)
          {
            int index = (i*(NUM_VERTS_Y-1) + j)*6;
            gIndices[index++] = j*NUM_VERTS_X+i;
            gIndices[index++] = j*NUM_VERTS_X+i+1;
            gIndices[index++] = (j+1)*NUM_VERTS_X+i+1;

            gIndices[index++] = j*NUM_VERTS_X+i;
            gIndices[index++] = (j+1)*NUM_VERTS_X+i+1;
            gIndices[index++] = (j+1)*NUM_VERTS_X+i;
            //cout<<"Indices: "<<j*NUM_VERTS_X+i<<"\t"<<j*NUM_VERTS_X+i+1<<(j+1)*NUM_VERTS_X+i+1<<"\t"<<j*NUM_VERTS_X+i<<"\t"<<(j+1)*NUM_VERTS_X+i+1<<"\t"<<(j+1)*NUM_VERTS_X+i<<endl;
          }
        }
        return CreateMeshFromData(m_vertices, gIndices, totalTriangles, totalVerts);
      }

			//Simulates the time dt in seconds with the specified number of substeps
			void Step(double dt, int substeps)
			{
				if(!m_dynamicsWorld)
        {
          cerr<<"m_dynamicsWorld not defined"<<endl;
          return;
        }
				double fixedstepsize = dt/substeps;
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
#endif
