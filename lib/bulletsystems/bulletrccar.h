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
   * A simple rear-drive Rccar model derived from Vehicle Demo using Bullet Physics Engine
   *
	 * This car model uses the full 6DOF state for simulation. To keep compatibility with existing 2D rccar model, 
	 * this class only provides the state output as \f$ \bm x = (x,y,\theta, v) \in\mathbb{R}^4\f$.
	 * The controls for the car are given as desired steering angle and desired body velocity. 
	 * These controls are converted into true controls for the bullet rccar i.e driving torque and actual steering angle 
	 * through proportional control. The gains for the controller are input as parameters to the system
	 * #TODO Use full state of the car
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   * Modified: Gowtham Garimella ggarime(at)jhu.edu
   */
  class Bulletrccar : public Rccar 
  {
    public:

			/** Full car state
			 */
    struct CarState{
      btTransform cartransform;
      btVector3 carlinearvel;
      btVector3 carangularvel;
    };

		/** Constructor. Each bullet system takes in a world in which it lies. 
		 * @param m_world 		World class in which the current system lies. Only systems in one world can interact with each other
		 * @param zs_					Vector of car heights along the trajectory used for display purposes.
		*/
    Bulletrccar(BulletWorld& m_world, vector<double> *zs_ = 0);

		/** Destructor: Removes the car and its accessories from the world and deletes them
		*/
    ~Bulletrccar()
    {
      if(m_carChassis) 
      {
        //m_world.m_dynamicsWorld->removeRigidBody(m_carChassis);

        delete m_carChassis;
      }

      delete m_vehicleRayCaster;

      m_world.m_dynamicsWorld->removeVehicle(m_vehicle);
      delete m_vehicle;

      delete m_wheelShape;
      
      delete initialstate;

      
      //m_world.Reset();
    }
    
		/** Reimplementation of the step function specific to bullet rccar
		 */
    double Step(Vector4d &xb, double t, const Vector4d &xa,
                const Vector2d &u, double h, const VectorXd *p,
                Matrix4d *A = 0, Matrix42d *B = 0, Matrix4pd *C = 0);

    double Step_internalinput(Vector4d &xb, const Vector2d &u, 
                double h, const VectorXd *p = 0,
                Matrix4d *A = 0, Matrix42d *B = 0, Matrix4pd *C = 0);

    double Step_noise(Vector4d &xb, const Vector2d &u,
                      const Vector4d &w, double h,
                      const VectorXd *p = 0,Matrix4d *A = 0, Matrix42d *B = 0, Matrix4pd *C = 0, Matrix4d *D = 0);

    bool reset(const Vector4d &x, double t = 0);

		/** Set the intial state of the car to given full state. This overrides the reset state of the car and allows to use full state
		*  #TODO Use full stat of the car directly
		*/
    void setinitialstate(const CarState &inputstate, Vector4d &x);
    void setinitialstate(Vector4d &x);///< Set the car initial state as the current state

		/** Reimplementation of Noise matrix for specific system
		 */
    bool NoiseMatrix(Matrix4d &Q, double t, const Vector4d &x, const Vector2d &u, double dt, const VectorXd *p);

    int rightIndex;///< right of the car index in terms of x,y,z assigns the coordinate sys of car
    int upIndex;
    int forwardIndex;

    double carmass;///< Mass of the car
    btVector3 car_halfdims;///< half dimensions of car width, height, length

    double	gEngineForce ;///< Engine Torque
    double	gBreakingForce;///< Breaking Force on the vehicle

    double maxEngineForce;///< Maximum Engine force that can be applied
    double	maxBreakingForce;///< Not Used 

    double	gVehicleSteering;///< Vehicle Steering angle
    double	steeringClamp;///< Clamp on Steering angle
    double velocityClamp;///< Clamp on velocity of car
    double	wheelRadius;///< Radius of wheel
    double	wheelWidth;///< Width of wheel
    double	wheelFriction;///< wheel friction usually large value
    double	suspensionStiffness;///< Stiffness of suspension. If not enough, the car will sink into ground
    double	suspensionDamping;///< Damping on suspension
    double	suspensionCompression;///< Compression factor decides by how much the car will compress wrto external loads
    double	rollInfluence;//< Decides whether the car will topple or not
    btScalar suspensionRestLength;///< Rest length of suspension
    btScalar m_defaultContactProcessingThreshold;///< if contact goes above this value, it will process

    btVector3 wheelDirectionCS0;///<direction from the car towards the wheel contact point
    btVector3 wheelAxleCS;///<Wheel Axle Direction

    //Parameters for the car:
    double gain_cmdvelocity;///<Gain parameter on the commanded velocity. Can be used to convert radio commands into true velocity commands
    double kp_torque;///<Proportional gain parameter for commanding torque to get to a desired velocity
    double kp_steer;///<Steering proportional gain. Should be between 0 and 1
    double initialz;///< Initial height of the car.

    btTransform offsettrans;///< To Account for the difference to coordinate system Usage: worldpose_inregcoordsys = offsettrans.inv()*worldpose_bullet*offsettrans
    btTransform offsettransinv;///< The inverse of the offset transformation above

    //Additional Hidden state for publishing/rendering trajectories:
    vector<double> *zs;///< Height of the rccar only used for visualization only x and y are used in optimization
    int count_zs;///< Internal variable to find the zs count to fill the height

    //Hacky way of ensuring the vehicle_velocity is reset after every reset
    bool reset_drivevel;///< Used to make sure the velocity of the chassis is used just after reset instead of vehicle vel

    //Bullet classes for holding car
    btRigidBody* m_carChassis;///< Pointer to the car rigidbody
    btCollisionShape* chassisShape;///< Pointer to the car collision shape
    btRaycastVehicle*	m_vehicle;///< Pointer to the raycast vehicle (Does the constraint enforcing to convert the rigidbody to a vehicle)
    btCollisionShape*	m_wheelShape;///< Shape of the wheels used only for display purposes
    btRaycastVehicle::btVehicleTuning	m_tuning;///< Tuning parameters for the vehicle used only internally
    btVehicleRaycaster*	m_vehicleRayCaster;///< Raycaster which gives the point of contact of the vehicle wheels
    BulletWorld &m_world;//< Parent world to which the car belongs
    CarState *initialstate;///< Initialstate if passed is stored here
  };
}


#endif
