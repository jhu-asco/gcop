#ifndef GCOP_INS_H
#define GCOP_INS_H

#include "system.h"
#include "insmanifold.h"
#include "so3.h"
#include <limits>
#include <iostream>
#include <utility>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  typedef Matrix<double, 15, 1> Vector15d;
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 15, 6> Matrix15x6d;
  typedef Matrix<double, 15, 15> Matrix15d;
  typedef Matrix<double, 15, Dynamic> Matrix15Xd;

  typedef Matrix<double, 6, 3> Matrix63d;
  typedef Matrix<double, 6, 6> Matrix6d;
  typedef Matrix<double, 6, 15> Matrix6x15d;
  typedef Matrix<double, 6, Dynamic> Matrix6Xd;
  
  /**
   * An INS system
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Ins : public System<InsState, 15, 6> {        
  public:
  
  Ins();    
  
  virtual ~Ins();
  
  double Step(InsState &xb, double t, const InsState &xa, 
              const Vector6d &u, double h, const VectorXd *p = 0,
              Matrix15d *A = 0, Matrix15x6d *B = 0, Matrix15Xd *C = 0);
  
  bool Noise(Matrix15d &Q, double t, const InsState &x, const Vector3d &u, 
             double dt, const VectorXd *p = 0);
  
  double sv; ///< gyro bias white noise stdev (spectral density)
  double su; ///< gyro bias rate-of-change white noise stdev (spectral density)
  double sa; ///< acceleration bias rate-of-change white noise stdev (spectral density)  

  double sra; ///< acceleration measurement noise

  Vector3d g0;   ///< gravity vector
  };  
}

#endif
