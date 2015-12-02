#ifndef GCOP_DISK_H
#define GCOP_DISK_H

#include <Eigen/Dense>

namespace gcop {
  
  class Disk {
  public:
    Disk(const Eigen::Vector2d &o, double r);
    
    Eigen::Vector2d o;  ///< origin
    double r;           ///< radius  
  };
}



#endif
