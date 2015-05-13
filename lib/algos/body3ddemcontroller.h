#ifndef GCOP_BODY3DDEMCONTROLLER_H
#define GCOP_BODY3DDEMCONTROLLER_H


#include "dsl/gridsearch.h"
#include "dsl/distedgecost.h"
#include "body3davoidcontroller.h"
#include "pqpdem.h"
#include "controller.h"

namespace gcop {

  using namespace dsl;  
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 5, 1> Vector5d;
   
  typedef PqpDem<Body3dState, 12, 6> Body3dPqpDem;
  
  /**
   * Controller combining rigid body stablization and gyroscopic avoidance
   * This basically combines two existing controller from DGC.
   * The controller has 5 parameters s containing: rot_Kp, rot_Kd, trans_Kp, trans_Kd, obst_K
   * 
   * Author: Marin Kobilarov, 2015
   */
  template <int nu = 6> class Body3dDemController : public Controller<Body3dState, Matrix<double, nu, 1>, Vector5d, Body3dState > {
  public:
  
  typedef Matrix<double, nu, 1> Vectorcd;
  
  /**
   * BODY3D PD controller
   *
   * @param sys system
   * @param xd desired state (optional, set to origin by default)
   * @param ad desired acceleration (optional, set to zero by default)
   * @param con obstacle constraint (optional)
   */
  Body3dDemController(const Body3d<nu> &sys,
                      const Body3dState &x0,
                      Body3dState *xf,
                      Body3dPqpDem &dem,
                      double vd = 20);
  
  virtual bool Set(Vectorcd &u, double t, const Body3dState &x);
  
  virtual bool SetParams(const Vector5d &s);
  
  virtual bool SetContext(const Body3dState &x);
  
  //  virtual bool SetContext(const Body3dPqpDem &c);
  
  void GetTraj(vector<Body3dState> &xds, const Dem &dem, GridSearch &gdsl, const Body3dState &x0, const Body3dState &xf, double vd);
  
  const Body3dState &x0;               ///< start state
  Body3dPqpDem &dem;                   ///< gidital elevation map
  Body3dAvoidController<nu> localCtrl; ///< local avoidance controller
  GridSearch dsl;                      ///< global grid search
  vector<Body3dState> xds;             ///< computed sequence of waypoints
  double vd;                           ///< desired forward velocity
  int j;                               ///< waypoint index
  double wpRadius;                     ///< waypoint switch radius
  };

  
  template <int nu>   
    void Body3dDemController<nu>::GetTraj(vector<Body3dState> &xds, const Dem &dem, GridSearch &gdsl, const Body3dState &x0, const Body3dState &xf, double vd) {
    
    int i0,j0,ig,jg;
    dem.Point2Index(i0, j0, x0.second[0], x0.second[1]);
    dem.Point2Index(ig, jg, xf.second[0], xf.second[1]);
    gdsl.SetStart(j0, i0);
    gdsl.SetGoal(jg, ig);
    GridPath path, optPath;
    gdsl.Plan(path);
    gdsl.OptPath(path, optPath, 2);
    for (int i = 0; i < optPath.count; ++i) {
      Body3dState x = xf;
      dem.Index2Point(x.second[0], x.second[1], optPath.pos[2*i+1], optPath.pos[2*i]);
      x.second[2] = x0.second[2];
      
      // if not last point
      Vector3d v(0,0,0);
      if (i > 0) {
        Vector3d pa;
        Vector3d pb;
        dem.Index2Point(pa[0], pa[1], optPath.pos[2*i-1], optPath.pos[2*i-2]);
        dem.Index2Point(pb[0], pb[1], optPath.pos[2*i+1], optPath.pos[2*i]);
        pa[2] = x0.second[2];
        pb[2] = x0.second[2];
        v = pb - pa;
        v = v/v.norm();
        v = v*vd;
      }
      x.second.tail<3>() = v;      
      xds.push_back(x);
    }
  }


  template <int nu>   
    Body3dDemController<nu>::Body3dDemController(const Body3d<nu> &sys,
                                                 const Body3dState &x0,
                                                 Body3dState *xf,
                                                 Body3dPqpDem &dem,
                                                 double vd) :
    x0(x0), dem(dem), localCtrl(sys, xf, 0, &dem), dsl(dem.dem.nj, dem.dem.ni, new DistEdgeCost(), dem.dem.data), vd(vd), j(0), wpRadius(10)
    {
    }

  template <int nu>   
    bool Body3dDemController<nu>::Body3dDemController::Set(Vectorcd &u, 
                                                           double t, 
                                                           const Body3dState &x)
    {    
      if (!xds.size()) {
        cout << "[W] Body3dDemController::Set: no desired states set yet!" << endl;
        return false;
      }

      assert(j < xds.size());
      
      Vector3d d = x.second.head<3>() - xds[j].second.head<3>();//localCtrl.stabCtrl.xd->second.head(3);
      if (d.norm() < wpRadius && j < xds.size() - 1)
        ++j;
      
      // set current waypoint
      localCtrl.stabCtrl.xd = &xds[j];

      // get local control
      localCtrl.Set(u, t, x);
      this->localCtrl.stabCtrl.sys.U.Clip(u);
                  
      return true;
    }
  
  template <int nu>   
    bool Body3dDemController<nu>::SetParams(const Vector5d &s) {
    return localCtrl.SetParams(s);
  }
  
  template <int nu>   
    bool Body3dDemController<nu>::SetContext(const Body3dState &xf) {

    // make sure it is free
    Matrix<double, 1, 1> g;
    dem(g, 0, xf);
    if (g[0] > 0)
      return false;

    xds.clear();
    // compute a new sequence of waypoint
    GetTraj(xds, dem.dem, dsl, x0, xf, vd);
    // reset index
    j = 0;

    return true;
  }

};

#endif
