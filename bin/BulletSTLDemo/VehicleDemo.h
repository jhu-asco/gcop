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
#ifndef STL_DEMO_H
#define STL_DEMO_H

#include "GlutDemoApplication.h"

#include <stdio.h>

#include <stdlib.h>

#include "utils.h"
#include "bulletworld.h"
#include "bulletrccar.h"

using namespace std;
using namespace gcop;
using namespace Eigen;

///VehicleDemo shows loading of a simple mesh into bullet along with a rccar:
class VehicleDemo : public GlutDemoApplication
{
  protected:
    Bulletrccar *brccar;
    BaseSystem *base_brccar;

	public:

	float		m_cameraHeight;

	float	m_minCameraDistance;
	float	m_maxCameraDistance;

  float gVehicleSteering;
  float gVehicleVelocity;
  float	steeringIncrement;//Steering Increment for the steering angle
  float velocityIncrement;//Velocity Increment for vehicle

  BulletWorld *world;

  VehicleDemo();

	virtual ~VehicleDemo();

	virtual void	clientResetScene();

  void clientMoveAndDisplay(void);

	virtual void displayCallback();

	virtual void updateCamera();

	virtual void specialKeyboard(int key, int x, int y);

	virtual void specialKeyboardUp(int key, int x, int y);

	virtual void keyboardCallback(unsigned char key, int x, int y);
	
	void renderme();

	void initPhysics();

	static DemoApplication* Create()
	{
		VehicleDemo* demo = new VehicleDemo();
		demo->myinit();
		demo->initPhysics();
		return demo;
	}
};

#endif //VEHICLE_DEMO_H


