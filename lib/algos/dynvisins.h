#ifndef GCOP_DYNVISINS_H
#define GCOP_DYNVISINS_H

#include <cmath>
#include <cstdio>
#include <iostream>

#include <fstream>
#include "body3d.h"
#include <Eigen/Dense>

#include "ceres/ceres.h"
#include "ceres/rotation.h"

namespace gcop {

using namespace std;
using namespace Eigen;

/**
 * Dynamic Visual-inertial bundle adjustment. The system trajectory is composed of
 * N segments, where a "segment" is defined as the time interval b/n two
 * camera measurements. There are at least one IMU measurement within each segment.
 * System dynamics can also be enforced, with default implementation using a constant
 * velocity model. 
 *
 * The trajectory is C-2 smooth (twice differentiable) and parametrized using a cubic spline
 * on each segment (C-4 smooth), so that the velocities also match at the segment endpoinds. 
 * This parametrization is achieved using the sequence (r0,p0,dr0,v0,r1,p1,dr1,v1,...,rN,pN,drN,vN)
 * where r are exponential coordinates, p is the position, dr is the rate of change of r, and 
 * v is the rate of change of p (i.e. the spatial velocity). 
 *
 * Adding likelihoods for IMU data, dynamics, or a given prior are optional.
 * Camera feature data is required since it forms the basis for the representation and 
 * optimization. All likelihood factors are computed in closed form without approximations.
 * 
 * In the current implementation, biases are not optimized.
 * 
 * Author: Marin Kobilarov
 */
class DynVisIns {
public:
  
  ceres::Problem problem;

  vector<Vector3d> ls;     ///< 3d points
  vector<Body3dState> xs;  ///< sequence of states

  double *v;             ///< the full ceres optimization vector

  vector<Vector2d> zs;   ///< all perpsective measurements
  vector<Vector3d> lus;  ///< all unit-spherical measurements (points in the IMU body frame)
  vector<int> zInds;     ///< measurements points indices (into ls)
  vector<int> zCamInds;  ///< measurements camera indices (into xs)

  double pxStd;          ///< std dev of pixels
  double sphStd;         ///< induced std dev on spherical measurements
  
  double dwStd;           ///< standard deviation for white noise angular accel
  double dvStd;           ///< standard deviation for white noise transl accel
  
  double wStd;            ///< gyro measurement stdev (e.g. from gyro bias stdev)
  double aStd;            ///< acc measurement stdev

  Vector3d g0;            ///< gravity vector

  double fx, fy, cx, cy;   ///< the camera intrinsic parameters
  Matrix<double, 2, 3> K;  ///< camera instrinsic matrix, for convenience

  Matrix3d Ric;            ///< from IMU to Camera rotation

  double t;     ///< global time
  double tc;    ///< time at last camera frame

  Body3dState x0;  ///< ins sequence start state (used as a local reference frame)

  Vector3d bg;  ///< acceleration bias
  Vector3d ba;  ///< gyro bias

  bool useImu;     ///< process IMU?
  bool useCam;     ///< process cam? 
  bool useDyn;     ///< use dynamics?

  bool usePrior;   ///< whether to enforce prior using x0

  bool optBias;    ///< to optimize over biases?

  bool sphMeas;    ///< use spherical measurements instead of standard pixel perspective measurements?


  vector<double> dts;           ///< camera segment delta times
  vector<vector<double> > tss;  ///< local times (within each segment) at which IMU measurements arrived
  vector<vector<Vector3d> > wss;  ///< accumulated IMU gyro readings from last cam frame
  vector<vector<Vector3d> > ass;  ///< accumulated IMU acc readings from last cam frame
  
  DynVisIns();
  virtual ~DynVisIns();

  /**
   * Convert from stl/Eigen data structures to ceres optimization vector
   * @param v ceres optimization vector
   * @return true if success
   */
  bool ToVec(double *v) {
    Vector12d c;
    for (int i = 0; i < xs.size(); ++i) {
      FromState(c, xs[i]);
      memcpy(v + 12*i, c.data(), 12*sizeof(double));
    }

    int i0 = 12*xs.size(); 
    for (int i = 0; i < ls.size(); ++i) {
      memcpy(v + i0 + 3*i, ls[i].data(), 3*sizeof(double));
    }

    if (optBias) {
      memcpy(v + 12*xs.size() + 3*ls.size(), bg.data(), 3*sizeof(double));
      memcpy(v + 12*xs.size() + 3*ls.size() + 3, ba.data(), 3*sizeof(double));
    }    
    return true;
  }
  

  /**
   * Convert from ceres optimization vector to stl/Eigen data structures
   * @param v ceres optimization vector
   * @return true if success
   */
  bool FromVec(const double *v) {
    for (int i = 0; i < xs.size(); ++i) {
      ToState(xs[i], Vector12d(v + 12*i));
    }

    int i0 = 12*xs.size(); 
    for (int i = 0; i < ls.size(); ++i) {
      ls[i] = Vector3d(v + i0 + 3*i);
    }

    if (optBias) {
      bg = Vector3d(v + 12*xs.size() + 3*ls.size());
      ba = Vector3d(v + 12*xs.size() + 3*ls.size() + 3);
    }

    return true;
  }

  static void ToState(Body3dState &x, 
                      const Vector3d &r,
                      const Vector3d &p,
                      const Vector3d &dr,
                      const Vector3d &v) {
    SO3::Instance().exp(x.R, r);
    x.p = p;
    Matrix3d D;
    SO3::Instance().dexp(D, -r);
    x.w = D*dr;
    x.v = v;
  }

  static void ToState(Body3dState &x, const Vector12d &c) {
    ToState(x, c.head<3>(), c.segment<3>(3), c.segment<3>(6), c.tail<3>());
  }


    
  static void FromState(Vector3d &r, Vector3d &p, Vector3d &dr, Vector3d &v,
                        const Body3dState &x) {
    SO3::Instance().log(r, x.R);
    p = x.p;
    Matrix3d D;
    SO3::Instance().dexpinv(D, -r);
    dr = D*x.w;
    v = x.v;
  }

  static void FromState(Vector12d &c, const Body3dState &x) {
    //    FromState(c.head<3>(), c.segment<3>(3), c.segment<3>(6), c.tail<3>(), x);
    
    Vector3d r;
    SO3::Instance().log(r, x.R);
    c.head<3>() = r;
    c.segment<3>(3) = x.p;
    Matrix3d D;
    SO3::Instance().dexpinv(D, -r);
    c.segment<3>(6) = D*x.w;
    c.tail<3>() = x.v;
  }

  
  /**
   * Estimate trajectory and points
   * @return true on success
   */
  bool Compute();

  /**
   * Process camera features
   * @param t time
   * @param zs perspective features
   * @param zInds features indices
   * @return true on success
   */
  bool ProcessCam(double t, const vector<Vector2d> &zs, const vector<int> &zInds);

  /**
   * Process IMU measurement
   * @param t time
   * @param w gyro measurement
   * @param a acceleration measurement
   * @return true on success
   */
  bool ProcessImu(double t, const Vector3d &w, const Vector3d &a);

  /**
   * Generate synthetic cam and IMU data and store in a provided "true" system tvi
   * @param tvi the true VI system
   * @param ni number of IMU measurements per camera segment
   */
  bool GenData(DynVisIns &tvi, int ni = 2);

  bool LoadFile(const char* filename);
};

}

#endif
