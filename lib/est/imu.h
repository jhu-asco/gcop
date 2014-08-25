#ifndef GCOP_IMU_H
#define GCOP_IMU_H

#include "system.h"
#include "imumanifold.h"
#include "so3.h"
#include <limits>
#include <iostream>
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 9, 1> Vector9d;
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 9, 3> Matrix93d;
  typedef Matrix<double, 9, 6> Matrix96d;
  typedef Matrix<double, 9, 9> Matrix9d;
  typedef Matrix<double, 9, Dynamic> Matrix9Xd;

  typedef Matrix<double, 6, 3> Matrix63d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 6, 9> Matrix69d;
  typedef Matrix<double, 6, Dynamic> Matrix6Xd;
  
  /**
   * An IMU system
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Imu : public System<ImuState, 9, 3> {        
  public:
  
  Imu();
  
  virtual ~Imu();
  
  double Step(ImuState &xb, double t, const ImuState &xa, 
              const Vector3d &u, double h, const VectorXd *p = 0,
              Matrix9d *A = 0, Matrix93d *B = 0, Matrix9Xd *C = 0);
  
  bool Noise(Matrix9d &Q, double t, const ImuState &x, const Vector3d &u, 
             double dt, const VectorXd *p = 0);
  
  double sv; ///< gyro bias white noise stdev (spectral density)
  double su; ///< gyro bias rate-of-change white noise stdev (spectral density)
  double sa; ///< acceleration bias rate-of-change white noise stdev (spectral density)

  };  
}

#endif
