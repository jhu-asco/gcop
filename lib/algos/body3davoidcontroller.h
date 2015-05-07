#ifndef GCOP_BODY3DAVOIDCONTROLLER_H
#define GCOP_BODY3DAVOIDCONTROLLER_H

#include "body3dcontroller.h"
#include "gavoidcontroller.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  /**
   * Controller combining rigid body stablization and gyroscopic avoidance
   * This basically combines two existing controller from DGC.
   * The controller has 5 parameters s containing: rot_Kp, rot_Kd, trans_Kp, trans_Kd, obst_K
   * 
   * Author: Marin Kobilarov, 2015
   */
  template <int nu = 6> 
    class Body3dAvoidController : public Controller<Body3dState, Matrix<double, nu, 1>, Matrix<double, 5, 1>, Body3dState > {
  public:
  typedef Matrix<double, nu, 1> Vectorcd;
  typedef Matrix<double, 5, 1> Vector5d;
  
  /**
   * BODY3D PD controller
   *
   * @param sys system
   * @param xd desired state (optional, set to origin by default)
   * @param ad desired acceleration (optional, set to zero by default)
   * @param con obstacle constraint (optional)
   */
  Body3dAvoidController(const Body3d<nu> &sys,
                        Body3dState *xd = 0,
                        Vector6d *ad = 0,
                        Body3dConstraint *con = 0);
  
  virtual bool Set(Vectorcd &u, double t, const Body3dState &x);

  virtual bool SetParams(const Vector5d &s);

  virtual bool SetContext(const Body3dState &c);

  Body3dController<nu> stabCtrl;
  GavoidController *avoidCtrl;
  
  };
  
  template <int nu>   
  Body3dAvoidController<nu>::Body3dAvoidController(const Body3d<nu> &sys,
                                                  Body3dState *xd, 
                                                  Vector6d *ad,
                                                  Body3dConstraint *con) :
    stabCtrl(sys, xd, ad)
    {
      if (con) 
        avoidCtrl = new GavoidController(*con);
      else
        avoidCtrl = 0;
    }
  
  template <int nu>   
    bool Body3dAvoidController<nu>::Body3dAvoidController::Set(Vectorcd &u, 
                                                               double t, 
                                                               const Body3dState &x)
    {
      stabCtrl.Set(u, t, x);      
      if (avoidCtrl) {
        Vector6d uo;
        avoidCtrl->Set(uo, t, x);
        u += uo; // add translational forces for obstacle avoidance
      }
      return true;
    }
  
  template <int nu>   
    bool Body3dAvoidController<nu>::SetParams(const Vector5d &s) {
    stabCtrl.Kp[0] = s[0];
    stabCtrl.Kp[1] = s[0];
    stabCtrl.Kp[2] = s[0];
    stabCtrl.Kd[0] = s[1];
    stabCtrl.Kd[1] = s[1];
    stabCtrl.Kd[2] = s[1];
    stabCtrl.Kp[3] = s[2];
    stabCtrl.Kp[4] = s[2];
    stabCtrl.Kp[5] = s[2];
    stabCtrl.Kd[3] = s[3];
    stabCtrl.Kd[4] = s[3];
    stabCtrl.Kd[5] = s[3];
    
    if (avoidCtrl)
      avoidCtrl->k = s[4];

    return true;
  }

  template <int nu>   
    bool Body3dAvoidController<nu>::SetContext(const Body3dState &c) {
    
    if (avoidCtrl) { // check for feasibility
      Matrix<double, 1, 1> g;
      (avoidCtrl->con)(g, 0, c);

      //      cout << "CHecking" << c.second.transpose() << " g=" << g << endl;
      
      if (g[0] > 0)
        return false;
    }
        
    stabCtrl.xd = (Body3dState*)&c;

    return true;
  }

};

#endif
