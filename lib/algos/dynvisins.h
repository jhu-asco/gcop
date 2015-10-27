#ifndef GCOP_DYNVISINS_H
#define GCOP_DYNVISINS_H

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>


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
  
  struct Point {
    Vector3d l;              ///< the point 3d coordinates
    Matrix3d P;              ///< the covariance matrix of the point 3d coordinates
    bool usePrior;           ///< should a prior be set based on the point covariance?
    map<int, Vector2d> zs;   ///< map of feature measurements indexed by camera id
  };

  struct Camera {
    Body3dState x;              ///< the point 3d coordinates
    vector<int> pntIds;         ///< list of points observed by this camera

    double dt;                 ///< delta t to next camera

    vector<double> ts;        ///< local times (within each segment) at which IMU measurements arrived
    vector<Vector3d> ws;      ///< accumulated IMU gyro readings from last cam frame
    vector<Vector3d> as;      ///< accumulated IMU acc readings from last cam frame    
  };
  
  map<int, Point> pnts;   ///< all points
  
  int camId;              ///< current camera id (incremented after a frame is added)
  int camId0;             ///< starting camera id (pointing to begining of window)
  map<int, Camera> cams;  ///< cameras

  int maxCams;            ///< max length of camera sequence (0 by default indicating no limit)

  ceres::Problem* problem;  ///< the ceres problem
  bool ceresActive;       ///< whether ceres has been called on the current opt vector

  double *v;             ///< the full ceres optimization vector
  std::vector<int> *l_opti_map;             ///< maps pnts in optimization vector back to pntIds
  int num_opti_cams;     ///< number of cameras used in the last optimization
  int n_good_pnts;        ///< number of good points in the last optimization

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

  Body3dState x0;  ///< ins sequence start state

  Vector3d bg;  ///< acceleration bias (assumed fixed during optimization window)
  Vector3d ba;  ///< gyro bias (assumed fixed during optimization window)

  bool useImu;     ///< process IMU?
  bool useCam;     ///< process cam? 
  bool useDyn;     ///< use dynamics?
  bool useAnalyticJacs;     ///< use analytic jacobians?
  bool useCay;      ///< use cayley map instead of exponential map?

  bool usePrior;   ///< whether to enforce prior using x0
  bool useFeatPrior;   ///< whether to enforce feature prior

  bool optBias;    ///< to optimize over biases?

  bool sphMeas;    ///< use spherical measurements instead of standard pixel perspective measurements? (false by default)
  double maxIterations; ///< maximum number of solver iterations
   
  DynVisIns();
  virtual ~DynVisIns();

  
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
   * Remove a camera frame from sequence
   * @param id cam id
   * @param pnt_zs_removed optional output argument which is holds ids of points which had measurements removed yet are still present in the optimization.
   * @return true on success
   */
  bool RemoveCamera(int id, std::set<int>* pnt_zs_removed = NULL);

  /**
   * Remove a point
   * @param id point id
   * @return true on success
   */
  bool RemovePoint(int id);

  /**
   * Reset prior on initial state x0 and its covariance to the 
   * mean and covariance of camera#id  in the current sequence. 
   * This is called internally by ProcessCam in case when the oldest
   * frame needs to be dropped to stay within maxCams window; the prior is then
   * reset to the second camera, before removing the first.
   * @param id camera id to which the prior will be reset
   * @param pnts point ids on which the prior will be reset
   * @return true on success
   */
  bool ResetPrior(int id, std::set<int>* pnts = NULL);

  /**
   * Process IMU measurement
   * @param t time
   * @param w gyro measurement
   * @param a acceleration measurement
   * @return true on success
   */
  bool ProcessImu(double t, const Vector3d &w, const Vector3d &a);


  bool MakeFeature(Vector2d &z, const Body3dState &x, const Vector3d &l);

  bool MakeFeatures(vector<Vector2d> &zs, 
                    vector<int> &pntIds,
                    const Body3dState &x, 
                    const map<int, Point> &pnts);

  /**
   * Generate synthetic cam and IMU data and store in a provided "true" system tvi
   * @param tvi the true VI system
   * @param ns number segments
   * @param np number of points
   * @param ni number of IMU measurements per camera segment
   */
  bool GenData(DynVisIns &tvi, int ns, int np, int ni = 2);


  bool SimData(DynVisIns &tvi, int ns, int np, int ni = 2, double dt = .1);

  bool LoadFile(const char* filename);

  /**
   * Convert from stl/Eigen data structures to ceres optimization vector
   * @param v ceres optimization vector
   * @return true if success
   */
  bool ToVec(double *v, vector<int>* l_map) {
    Vector12d c;
    map<int, Camera>::iterator camIter;
    int i = 0;
    for (camIter = cams.begin(); camIter != cams.end(); ++camIter, ++i) {
      FromState(c, camIter->second.x);
      memcpy(v + 12*i, c.data(), 12*sizeof(double));
    }

    int i0 = 12*cams.size(); 
    map<int, Point>::iterator pntIter;
    i=0;
    for (pntIter = pnts.begin(); pntIter != pnts.end(); ++pntIter) {
      // Only consider points with at least two observations
      if(pntIter->second.zs.size() > 1)
      {
        memcpy(v + i0 + 3*i, pntIter->second.l.data(), 3*sizeof(double));
        l_map->at(i) = pntIter->first;
        ++i;
      }
    }

    //    if (optBias) {
    //      memcpy(v + 12*xs.size() + 3*ls.size(), bg.data(), 3*sizeof(double));
    //      memcpy(v + 12*xs.size() + 3*ls.size() + 3, ba.data(), 3*sizeof(double));
    //    }    
    return true;
  }
  

  /**
   * Convert from ceres optimization vector to stl/Eigen data structures
   * @param v ceres optimization vector
   * @return true if success
   */
  bool FromVec(const double *v, const vector<int>* l_map) {
    //    xs.resize(cams.size();
    map<int, Camera>::iterator camIter;
    int i=0;
    for (camIter = cams.begin(); camIter != cams.end(); ++camIter, ++i) {
      ToState(camIter->second.x, Vector12d(v + 12*i));
      //      xs[i] = camIter->second.x;
    }

    int i0 = 12*cams.size(); 
    //map<int, Point>::iterator pntIter;
    //i = 0;
    //for (pntIter = pnts.begin(); pntIter != pnts.end(); ++pntIter, ++i) {
    for (int i = 0; i < l_map->size(); ++i) {
      pnts[l_map->at(i)].l = Vector3d(v + i0 + 3*i);
    }

    //    if (optBias) {
    //      bg = Vector3d(v + 12*xs.size() + 3*ls.size());
    //      ba = Vector3d(v + 12*xs.size() + 3*ls.size() + 3);
    //    }

    return true;
  }

  void ToState(Body3dState &x, 
                      const Vector3d &r,
                      const Vector3d &p,
                      const Vector3d &dr,
                      const Vector3d &v) {
    Matrix3d D;
    if(useCay)
    {
      SO3::Instance().cay(x.R, r);
      SO3::Instance().dcay(D, -r);
    }
    else
    {
      SO3::Instance().exp(x.R, r);
      SO3::Instance().dexp(D, -r);
    }
    x.p = p;
    x.w = D*dr;
    x.v = v;
  }

  void ToState(Body3dState &x, const Vector12d &c) {
    ToState(x, c.head<3>(), c.segment<3>(3), c.segment<3>(6), c.tail<3>());
  }


    
  void FromState(Vector3d &r, Vector3d &p, Vector3d &dr, Vector3d &v,
                        const Body3dState &x) {
    Matrix3d D;
    if(useCay)
    {
      SO3::Instance().cayinv(r, x.R);
      SO3::Instance().dcayinv(D, -r);
    }
    else
    {
      SO3::Instance().log(r, x.R);
      SO3::Instance().dexpinv(D, -r);
    }
    p = x.p;
    dr = D*x.w;
    v = x.v;
  }

  void FromState(Vector12d &c, const Body3dState &x) {
    //    FromState(c.head<3>(), c.segment<3>(3), c.segment<3>(6), c.tail<3>(), x);
    
    Vector3d r;
    Matrix3d D;
    if(useCay)
    {
      SO3::Instance().cayinv(r, x.R);
      SO3::Instance().dcayinv(D, -r);
    }
    else
    {
      SO3::Instance().log(r, x.R);
      SO3::Instance().dexpinv(D, -r);
    }
    c.head<3>() = r;
    c.segment<3>(3) = x.p;
    c.segment<3>(6) = D*x.w;
    c.tail<3>() = x.v;
  }

};

}

#endif
