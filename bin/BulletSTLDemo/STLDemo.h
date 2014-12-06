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

using namespace std;
using namespace gcop;

///STLDemo shows loading of a simple mesh into bullet
class STLDemo : public GlutDemoApplication
{
	public:

	float		m_cameraHeight;

	float	m_minCameraDistance;
	float	m_maxCameraDistance;

  char *filename;

  BulletWorld *world;

  STLDemo();

	virtual ~STLDemo();

	//virtual void	clientResetScene();
  void clientMoveAndDisplay(){};

	virtual void displayCallback();
	
	//void renderme();

	void initPhysics();

	static DemoApplication* Create()
	{
		STLDemo* demo = new STLDemo();
		demo->myinit();
		demo->initPhysics();
		return demo;
	}
};

#endif //VEHICLE_DEMO_H


