#ifndef GCOP_SPHERE_H
#define GCOP_SPHERE_H

#include <Eigen/Dense>

namespace gcop {
  
  class Sphere {
  public:
    Sphere(const Eigen::Vector3d &o, double r);
    
    Eigen::Vector3d o;      ///< origin
    double r;        ///< radius  
  }; 
}
#endif
