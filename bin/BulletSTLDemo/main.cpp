
#include "STLDemo.h"
#include "GlutStuff.h"
#include "GLDebugDrawer.h"
#include "btBulletDynamicsCommon.h"
GLDebugDrawer	gDebugDrawer;

int main(int argc,char** argv)
{

        STLDemo* stlDemo = new STLDemo;
        if(argc > 1)
          stlDemo->filename = argv[1];

        stlDemo->initPhysics(); 
		stlDemo->getDynamicsWorld()->setDebugDrawer(&gDebugDrawer);

        return glutmain(argc, argv,640,480,"Bullet STL Demo.", stlDemo);
}

