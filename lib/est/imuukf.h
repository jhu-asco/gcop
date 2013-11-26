#ifndef GCOP_IMUUKF_H
#define GCOP_IMUUKF_H

#include "ukf.h"
#include "imumodel.h"

namespace gcop {
  class ImuUKF : public UKF {
  public:
    ImuUKF(ImuModel &im);

    virtual ~ImuUKF();

    void Reset();
    
    bool Process(double t,
                 const Eigen::Vector3d &w,
                 const Eigen::Vector3d &a,
                 const Eigen::Vector3d &m);      
    
    double t;      ///< current time 

    Eigen::Vector3d u;   ///< temp storage for inputs
    Eigen::VectorXd z;   ///< temp storage for measurements    

  protected:
    bool updated;

  };
}

#endif
