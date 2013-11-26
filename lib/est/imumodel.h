#ifndef GCOP_IMUMODEL_H
#define GCOP_IMUMODEL_H

#include "model.h"

namespace gcop {

  /**
   * General IMU Model for attitude gcopimation using
   * gyro, magnetic, and accelerometer data.
   * 
   * Author: Marin Kobilarov -- Copyright (C) 2005
   */
  class ImuModel : public Model {
  public:
    /**
     * IMU Model
     * @param acc use accelerometer (true by default)
     * @param gmr use a global magnetic reference (defined by m0). Setting gmr=false is experimental 
     */
    ImuModel(bool acc = true,
             bool gmr = true);

    virtual ~ImuModel();
    
    virtual bool f(Eigen::VectorXd &xn,
                   const Eigen::VectorXd &x,
                   const Eigen::VectorXd &u);
    
    virtual bool h(Eigen::VectorXd &z, 
                   const Eigen::VectorXd &x);
    
    
    void I(double *y,
           const double *z,
           const double *x);
    
    /**
     * Checks if the inputs are within bounds
     * @param w gyro rates
     * @param a accelerations
     * @param m magnetometer readings
     * @return true if valid
     */
    bool IsValid(const Eigen::Vector3d &w,
                 const Eigen::Vector3d &a,
                 const Eigen::Vector3d &m) const;


    bool acc;        ///< use accelerometer
    bool gmr;        ///< use global magnetic reference defined by m0 (otherwise the reference would be the first measurement obtained when the filter starts)
    
    Eigen::Vector3d wn;     ///< gyro angular acceleration white noise
    Eigen::Vector3d wbn;    ///< gyro bias noise
    Eigen::Vector3d wbn0;   ///< initial gyro bias noise
    Eigen::Vector3d vn;     ///< velocity acceleration white noise 
    Eigen::Vector3d abn;    ///< accelerometer bias noise
    Eigen::Vector3d abn0;   ///< initial accelerometer bias noise

    Eigen::Vector3d mon;    ///< magnetometer observation noise

    Eigen::Vector3d g;      ///< gravity vector (default is (0,0,-9.81))
    Eigen::Vector3d m0;     ///< reference magnetic vector

    double dt;        ///< time-step ( set internally by Update)
  };
}

#endif
