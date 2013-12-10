#include "joint.h"
#include "se3.h"

using namespace gcop;

Joint::Joint() : a(Vector6d::Zero()), 
                 gp(Matrix4d::Identity()),
                 gc(Matrix4d::Identity()),
                 damping(0),
                 friction(0)
{
  Init();
}

Joint::Joint(const Vector6d &a,
             const Matrix4d &gp,
             const Matrix4d &gc,
             double damping,
             double friction) : 
  a(a), gp(gp), gc(gc), damping(damping), friction(friction)
{
  Init();
}
 
Joint::~Joint()
{
}

void Joint::Init()
{
  SE3::Instance().Ad(Ac, gc);
  SE3::Instance().inv(gpi, gp);
  SE3::Instance().inv(gci, gc);  
  S = Ac*a;
}

