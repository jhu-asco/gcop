// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_MBS_H
#define GCOP_MBS_H

#include "body3d.h"
#include "joint.h"
#include "mbsmanifold.h"
#include "function.h"
#include "mbscspace.h"
#include "mbstspace.h"
#include <exception>      // std::exception
#include <stdexcept>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  /**
   * General multi-body dynamics time-stepping. The discrete mechanics is
   * based on an implicit symplectic variational integrator which is second-order accurate
   * and momentum-balance preserving. 
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Mbs : public System<MbsState> {
    
  public:

    Mbs(int nb, int c, bool fixed = false, int np=0);

    virtual ~Mbs();

    void Init();
    
    double Step(MbsState& xb, double t, const MbsState& xa,
                const VectorXd &u, double h, const VectorXd *p = 0,
                MatrixXd *A = 0, MatrixXd *B = 0, MatrixXd *C = 0);


    double EulerStep(MbsState& xb, double t, const MbsState& xa,
                     const VectorXd &u, double h, const VectorXd *p = 0,
                     MatrixXd *A = 0, MatrixXd *B = 0);


    double HeunStep(MbsState& xb, double t, const MbsState& xa,
                    const VectorXd &u, double h, const VectorXd *p = 0,
                    MatrixXd *A = 0, MatrixXd *B = 0);

    double TrapStep(MbsState& xb, double t, const MbsState& xa,
                    const VectorXd &u, double h, const VectorXd *p = 0,
                    MatrixXd *A = 0, MatrixXd *B = 0);      


    void NE(VectorXd &e, const VectorXd &vdr, 
            MbsState &xb,               
            double t, const MbsState &xa, 
            const VectorXd &u, double h, const VectorXd *p = 0);


    /**
     * Inverse-dynamics function. Computes the residual of the
     * dynamical update from xa to xb given controls u. This is regarded as
     * the discrete-time equivalent of (b - B*u)
     * @param f residual
     * @param t time
     * @param x state
     * @param u controls
     */
    void ID(VectorXd &f,
            double t, const MbsState &x,
            const VectorXd &u);

    /**
     * Compute total bias b, i.e. such that M*a + b = B*u
     * @param b bias
     * @param t time
     * @param x state
     */
    void Bias(VectorXd &b,
              double t, const MbsState &x, const VectorXd *p = 0) const;

    /**
     * Discrete bias
     * @param b bias
     * @param t time
     * @param xb state
     * @param xa state
     * @param h time-step
     */
    void DBias(VectorXd &b,
                double t,
                const MbsState &xb, 
                const MbsState &xa, double h, const VectorXd *p = 0);    
    
    /**
     * Forward kinematics: given base pose and joint angles, update the poses of
     * all bodies in the graph.
     * @param x state
     */
    void FK(MbsState &x);

    /**
     * Kinematic-step. Given previous state xa, the next state xb is updated
     * by taking the joint angles from xa, and joint velocities from xb, updating the
     * joint angles in xb, and updating the xb body poses
     * @param xa previous state
     * @param xb next state
     * @param h time-step
     * @param impl whether to update configuration using the newly updated velocity
     */
    void KStep(MbsState &xb, const MbsState &xa, double h, bool impl = true);

    void NewtonEulerJacobian(MatrixXd &De, const MbsState &xb, const MbsState &xa, double h);  

    /**
     * Compute the mass matrix at a given state x
     * @param M mass matrix
     * @param x state
     */
    void Mass(MatrixXd &M, const MbsState &x) const;

    /**
     * Compute acceleration 
     * @param a acceleration
     * @param t time 
     * @param x state
     * @param u control inputs
     * @param h time-step
     */
    void Acc(VectorXd &a, double t, const MbsState& x, const VectorXd &u, double h, const VectorXd *p = 0);

    /**
     * Total resulting force on the system from external (e.g. gravity)
     * and internal (control) inputs
     * @param f total force f(x,u,t)
     * @param t time
     * @param x state
     * @param u controls
     * @param A jacobian Dxf
     * @param B jacobian Duf
     */
    virtual void Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u, 
                       MatrixXd *A = 0, MatrixXd *B = 0);
    
    void Rec(MbsState &x, double h);

    void ClampPose(MbsState &x, int i) const;

    void CheckLimits(MbsState &x, int i, double h) const;

    void GetImpulse(double f, const MbsState &x, int i, double h) const;

    void ClampVelocity(MbsState &x) const;

		void print(const MbsState &x) const; 

    int nb;                   ///< number of rigid bodies (including base body)

    bool fixed;               ///< whether base body is fixed

    vector<Body3d<>> links;   ///< links

    vector<Joint> joints;     ///< joints

    vector<Matrix6d> Ips;     ///< A'*I*A (nb) vector
    
    vector<int> pis;          ///< parent indexes
    vector<vector<int> > cs;  ///< child lists

    Vector3d ag;             ///< acceleration due to gravity (0, 0, -9.81) by default

    VectorXd damping;        ///< joint damping vector
    
    VectorXd lbK;            ///< joint lower bound K term
    VectorXd lbD;            ///< joint lower bound D term
    VectorXd ubK;            ///< joint upper bound K term
    VectorXd ubD;            ///< joint upper bound D term
    VectorXd fsl;             ///< spring force lower bound
    VectorXd fsu;             ///< spring force upper bound
    
    int basetype;            ///< type of base (used for interpreting base-body forces)
    static const int FIXEDBASE = 0;   ///< fixed base
    static const int FLOATBASE = 1;   ///< fully controllable floating base
    static const int AIRBASE = 2;     ///< specialized free-flying base with three torques and one lift force

    SE3 &se3;                 ///< singleton reference for performing SE(3) operations

    int method;                   ///< type of integration method employed (Euler by default)
    static const int EULER = 1;   ///< Euler time-stepping of dynamics
    static const int HEUN = 2;    ///< Heun's explicit 2nd order method
    static const int TRAP = 3;    ///< symplectic trapezoidal 2nd order method
    int iters;                    ///< max number of Newton iterations used in symplectic method

    string end_effector_name;///< Name of end effectors to which parameter forces(*p) are applied

    bool debug;
  };
}

#endif
