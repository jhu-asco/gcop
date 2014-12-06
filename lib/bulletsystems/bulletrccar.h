#ifndef GCOP_BULLETRCCAR_H
#define GCOP_BULLETRCCAR_H

#include "system.h"
#include "rccar.h"
#include <BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <btBulletDynamicsCommon.h>
#include "bulletworld.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 4, 2> Matrix42d;
  typedef Matrix<double, 4, Dynamic> Matrix4pd;
  //typedef Matrix<double, 6, 1> Vector6d;

   /**
   * A simple rear-drive Rccar model with 2nd order dynamics. 
   *
   * The state is
   * \f$ \bm x = (x,y,\theta, v) \in\mathbb{R}^4\f$ and controls are \f$ \bm u = (a, tan(phi) ) \in\mathbb{R}^5\f$ 
   * correspond to forward acceleration a and steering angle phi
   * This uses bullet to propagate the dynamics of the car. To be consistent with bullet 
   * #TODO Add steering angle to state since, the control has been changed to cmded steering angle instead of true steering angle
   * demos and OpenGL coordinate system, we use coordinate system of y being up, x being to right, z being to the back of the car. This may be weird but makes it simple for display purposes
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   * Modified: Gowtham Garimella
   */
  class Bulletrccar : public Rccar 
  {
  public:
    Bulletrccar(BulletWorld& m_world);

    ~Bulletrccar()
    {
      delete m_vehicleRayCaster;

      delete m_vehicle;

      delete m_wheelShape;
    }
    
    double Step1(Vector4d &xb, const Vector2d &u, 
                double h, const VectorXd *p,
                Matrix4d *A = 0, Matrix42d *B = 0, Matrix4pd *C = 0);

    double Step3(Vector4d &xb, const Vector2d &u,
                      const Vector4d &w, double h,
                      const VectorXd *p = 0,Matrix4d *A = 0, Matrix42d *B = 0, Matrix4pd *C = 0, Matrix4d *D = 0);

    bool reset(const Vector4d &x, double t = 0);

    bool NoiseMatrix(Matrix4d &Q, double t, const Vector4d &x, const Vector2d &u, double dt, const VectorXd *p);

    int rightIndex;///<right of the car index in terms of x,y,z assigns the coordinate sys of car
    int upIndex;
    int forwardIndex;

    float carmass;//Mass of the car
    btVector3 car_halfdims;//half dimensions of car width, height, length

    float	gEngineForce ;//Engine Torque
    float	gBreakingForce;//Breaking Force on the vehicle

    float	maxEngineForce;
    float	maxBreakingForce;

    float	gVehicleSteering;//Vehicle Steering angle
    float	steeringClamp;// Clamp on Steering angle
    float velocityClamp;//Clamp on velocity of car
    float	wheelRadius;//Radius of wheel
    float	wheelWidth;//Width of wheel
    float	wheelFriction;//wheel friction usually large value
    float	suspensionStiffness;//Stiffness of suspension. If not enough, the car will sink into ground
    float	suspensionDamping;//Damping on suspension
    float	suspensionCompression;//Compression factor decides by how much the car will compress wrto external loads
    float	rollInfluence;//Decides whether the car will topple or not
    btScalar suspensionRestLength;//Rest length of suspension
    btScalar m_defaultContactProcessingThreshold;// if contact goes above this value, it will process

    btVector3 wheelDirectionCS0;//direction from the car towards the wheel contact point
    btVector3 wheelAxleCS;//Wheel Axle Direction

    //Parameters for the car:
    float gain_cmdvelocity;
    float kp_torque;
    float kp_steer;//Ideally should be between 0 and 1(Can define a map to do this #TODO)
    float initialz;

    //Bullet classes for holding car
    btRigidBody* m_carChassis;
    btRaycastVehicle*	m_vehicle;
    btCollisionShape*	m_wheelShape;
    btRaycastVehicle::btVehicleTuning	m_tuning;
    btVehicleRaycaster*	m_vehicleRayCaster;
    BulletWorld &m_world;
    //btDynamicsWorld*		m_dynamicsWorld;
  };
}


#endif
