// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_SYSTEM_H
#define GCOP_SYSTEM_H

#include <Eigen/Dense>
#include <vector>
#include <assert.h>
#include "manifold.h"
#include <iostream>
#include "rn.h"

//#include "sensor.h"

/*! \mainpage Geometric Control, Optimization, and Planning (GCOP) software library
 *
 * \section Description
 * \subsection Intro
 *
 * This is a partial distribution focusing on optimal control functionality. 
 * Several example systems are included such as a simple car model, aerial vehicles including quadrotor and helicopter, 
 * and fixed-based and free-flying multi-body systems.
 *
 * In addition, a parameter-dependent discrete optimal control method is included for solving
 * sys id / adaptive control / parameter estimation problems. It is applied to a 
 * SLAM/bundle adjustment problem.
 * 
 * 
 * \subsection Installation
 * \subsubsection Build requirements
 *  g++, cmake, libeigen, opengl, glut
 *
 * \subsubsection Download
 * 
 * To obtain the source code, please email: marin(at)jhu.edu 
 *
 * \subsubsection Compilation
 * - To compile >: 
 * -- tar xfz gcop.tar.gz; cd gcop 
 * -- mkdir build; cd build; cmake ..; make
 * - To test an optimal control example>: cd bin; ./qrotortest, ./cartest, etc...
 * - To test a bundle adjustment example>: cd bin; ./batest
 *
 * \subsection API Reference
 * see docs/html
 *
 * \subsection Usage
 *
 *
 * \subsection Example
 *  see directory bin
 * \subsection Author
 * Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
 * \subsection Keywords
 * discrete optimal control, discrete mechanics, motion planning, nonlinear control
 */

namespace gcop {
  
  using namespace Eigen;
  
  
  /**
   * Control system interface for dynamic systems on manifolds. 
   * Defines the discrete dynamics and any structural/inertial system parameters 
   * such as body dimensions, mass matrices, etc... 
   *
   * Subclasses should provide implementation for the 
   * time-stepping Step function and optionally for process noise evolution
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic> class System {
    
  public:
  
  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;

  typedef Matrix<double, _nx, _nx> Matrixnd;
  typedef Matrix<double, _nx, _nu> Matrixncd;
  typedef Matrix<double, _nu, _nx> Matrixcnd;
  typedef Matrix<double, _nu, _nu> Matrixcd;
  
  //  typedef Matrix<double, Dynamic, 1> Vectormd;
  typedef Matrix<double, _np, _np> Matrixmd;
  typedef Matrix<double, _nx, _np> Matrixnmd;
  typedef Matrix<double, _np, _nx> Matrixmnd;
  
  /**
   * A system evolving on a state manifold manifold X with nu inputs and np parameters
   * @param X state manifold
   * @param nu number of control inputs
   * @param np number of parameters
   */
  System(Manifold<T, _nx> &X, int nu = 0, int np = 0);
  
  /**
   * Discrete dynamics update. The function computes the next system state xb given 
   * current state xa and applied forces u
   * @param xb resulting state
   * @param t current time
   * @param xa current state
   * @param u current control
   * @param h current time-step
   * @param p static parameters  (optional)
   * @param A jacobian w.r.t. x   (optional)
   * @param B jacobian w.r.t. u   (optional)
   * @param C jacobian w.r.t. p   (optional)
   */
  virtual double Step(T &xb, double t, const T &xa,
                      const Vectorcd &u, double h, const Vectormd *p = 0,
                      Matrixnd *A = 0, Matrixncd *B = 0, Matrixnmd *C = 0);

  /**
   * Discrete Dynamics Update. This function computes the next state xb using the internal input
   * state x_int and applied forces u
   * @param xb resulting state
   * @param t current time
   * @param u current control
   * @param h current time-step
   * @param p static parameters  (optional)
   * @param A jacobian w.r.t. x   (optional)
   * @param B jacobian w.r.t. u   (optional)
   * @param C jacobian w.r.t. p   (optional)
   */
  virtual double Step_internalinput(T &xb, const Vectorcd &u,
                       double h, const Vectormd *p = 0,
                       Matrixnd *A = 0, Matrixncd *B = 0, Matrixnmd *C = 0);

   /**
   * Discrete Dynamics update of internal state and no output
   * @param t current time
   * @param u current control
   * @param h current time-step
   * @param p static parameters  (optional)
   * @param A jacobian w.r.t. x   (optional)
   * @param B jacobian w.r.t. u   (optional)
   * @param C jacobian w.r.t. p   (optional)

   */
  virtual double Step_internaloutput(const Vectorcd &u,
      double h,const Vectormd *p = 0,
      Matrixnd *A = 0, Matrixncd *B = 0, Matrixnmd *C = 0);

  /** 
   * Additional step function with input noise of the same dimension as state x and propagates by one dynamic step
   * Noise in the system is modeled as x[i+1] = f(xi,ui,ti,p) + noisematrix(xi,ui,ti,p)*wi
   * The noise matrix is implemented as a virtual function and by default is just the identity
   * @param xb resulting state
   * @param u current control
   * @param h current time-step
   * @param w noise of same size as state
   * @param p static parameters  (optional)
   * @param A jacobian w.r.t. x   (optional)
   * @param B jacobian w.r.t. u   (optional)
   * @param C jacobian w.r.t. p   (optional)
   * @param D jacobian w.r.t  w   (optional)
   */
  virtual double Step_noise(T &xb, const Vectorcd &u,
                      const Vectornd &w, double h,
                      const Vectormd *p = 0,Matrixnd *A = 0, Matrixncd *B = 0, Matrixnmd *C = 0, Matrixnd *D = 0);

  /** Resets the internal state of the system to one specified
   * @param x State to reset to
   * @param t  time to reset to (optional)
   */
  virtual bool reset(const T &x, double t = 0);

  /**
   * Reconstruct the full state. Some states contain redundant parameters 
   * for efficiency and it is useful to reconstruct them to maintain a 
   * consisten state.
   * @param x state
   * @param h time-step
   */
  virtual void Rec(T &x, double h) {};

  /**
   * General print function to print the state
   */
  virtual void print(const T &x) const {};
  
  /**
   * For filtering purposes, it is often necessary to provide 
   * a discrete-time process noise covariance Q, typically used in
   * linear error propagation according to P = A*P*A' + Q
   * @param Q discrete-time process noise covariance
   * @param t start time
   * @param x start state
   * @param u control
   * @param dt time-step
   * @param p parameter (optional)
   * @return true if all arguments are feasible
   */
  virtual bool Noise(Matrixnd &Q, double t, const T &x, const Vectorcd &u, 
                     double dt, const Vectormd *p = 0);
  /**
   * Noise in the system is modeled as x[i+1] = f(xi,ui,ti,p) + noisematrix(xi,ui,ti,p)*wi
   * This function provides the necessary noise coefficient matrix
   * @param Q provides the noise matrix
   * @param t current time
   * @param x current state
   * @param u current control
   * @param h time-step
   * @param p parameter (optional)
   * @return true if all arguments are feasible
   */
  virtual bool NoiseMatrix(Matrixnd &Q, double t, const T &x, const Vectorcd &u, 
                     double h, const Vectormd *p = 0);


  /**
   * Convert a system state and controls to flat outputs for the system.
   * @param y the generated flat outputs
   * @param x given state
   * @param u given controls
   */
  virtual void StateAndControlsToFlat(VectorXd &y, const T &x, const Vectorcd &u);

  /**
   * Convert flat outputs and derivatives of flat outputs to the state and controls of the system.
   * @param y given flat outputs and derivatives of flat outputs
   * @param x generated state
   * @param u generated controls
   */
  virtual void FlatToStateAndControls(T &x, Vectorcd &u, const std::vector<VectorXd> &y);

  Manifold<T, _nx> &X;   ///< state manifold 
  Rn<_nu> U;             ///< control Euclidean manifold
  Rn<_np> P;             ///< parameter Euclidean manifold  
  T x;               ///< Internal State of the system
  int np;            ///< Number of Parameters
  double t;          ///< Internal time of the system
  
  };

  
  template <typename T, int _nx, int _nu, int _np> 
    System<T, _nx, _nu, _np>::System(Manifold<T, _nx> &X, int nu, int np) : 
    X(X), U(nu), P(np), np(np) {
  }

  /*
  template <typename T, typename Tu, int _nx, int _nu> 
    double System<T, Tu, _nx, _nu>::F(Vectornd &v, double t, const T& x,
                                    const Tu &u, double h,
                                     const VectorXd *p,
                                     Matrix<double, _nx, _nx> *A, Matrix<double, _nx, _nu> *B, 
                                     Matrix<double, _nx, Dynamic> *C) {
    std::cout << "[W] System::F: unimplemented!" << std::endl;
    return 0;
  }
  */

  template <typename T, int _nx, int _nu, int _np> 
    double System<T, _nx, _nu, _np>::Step(T& xb, double t, const T& xa,
                                          const Vectorcd &u, double h, const Vectormd *p,
                                          Matrixnd *A, Matrix<double, _nx, _nu> *B, 
                                          Matrix<double, _nx, _np> *C) {
    /*
    Vectornd v;
    if (_nx == Dynamic)
      v.resize(n);

    double d = F(v, t, xa, u, h, p, A, B, C);
    X.Retract(xb, xa, v); // in R^n this is just xb = xa + v
    Rec(xb, h);
    return d;
    */
    return 0;
  }
  template <typename T, int _nx, int _nu, int _np> 
    double System<T, _nx, _nu, _np>::Step_internalinput(T& xb, const Vectorcd &u, double h, const Vectormd *p,
                                          Matrixnd *A, Matrix<double, _nx, _nu> *B, 
                                          Matrix<double, _nx, _np> *C) {
    double result = this->Step(xb, this->t, this->x, u, h, p, A, B, C);
    this->x = xb;
    this->t += h;
    return result;
  }
  template <typename T, int _nx, int _nu, int _np> 
    double System<T, _nx, _nu, _np>::Step_internaloutput(const Vectorcd &u, double h, const Vectormd *p,
                                          Matrixnd *A, Matrix<double, _nx, _nu> *B, 
                                          Matrix<double, _nx, _np> *C) {
      T xb;
      double result = this->Step_internalinput(xb, u, h, p, A, B, C);
      return result;
  }

  template <typename T, int _nx, int _nu, int _np>
  double System<T, _nx, _nu, _np>::Step_noise(T &xb, const Vectorcd &u,
                      const Vectornd &w, double h, const Vectormd *p,
                      Matrixnd *A, Matrixncd *B, Matrixnmd *C, Matrixnd *D)
  {
    //double result = this->Step(xb, this->t, this->x, u, h, p, A, B, C);
    double result = this->Step_internalinput(xb, u, h, p, A, B, C);
    Matrixnd noise_matrix;
    this->NoiseMatrix(noise_matrix, this->t, this->x, u, h, p);
    this->X.Retract(xb, xb, (noise_matrix*w));
    this->x = xb;//Update internal state
    this->t = (this->t + h);
    return result;
  }

  template <typename T, int _nx, int _nu, int _np>
    bool System<T, _nx, _nu, _np>::reset(const T& x, double t)
    {
      this->x = x;//Copying state and time 
      this->t = t;
      return true;
    }


  template <typename T, int _nx, int _nu, int _np> 
    bool System<T, _nx, _nu, _np>::Noise(Matrixnd &Q, double t, const T &x, const Vectorcd &u, 
               double dt, const Vectormd *p) {
    std::cout << "[W] System::Noise: unimplemented! Subclasses should override." << std::endl;
    return false;
  }

  template <typename T, int _nx, int _nu, int _np> 
    bool System<T, _nx, _nu, _np>::NoiseMatrix(Matrixnd &Q, double t, const T &x, const Vectorcd &u, 
               double dt, const Vectormd *p) {
      Q.setIdentity();
    return true;
  }

  template <typename T, int _nx, int _nu, int _np> 
    void System<T, _nx, _nu, _np>::StateAndControlsToFlat(VectorXd &y, const T &x, 
               const Vectorcd &u) {
    y.resize(0);
    std::cout << "[W] System::StateAndControlsToFlat: unimplemented! Subclasses should override." << std::endl;
  }

  template <typename T, int _nx, int _nu, int _np> 
    void System<T, _nx, _nu, _np>::FlatToStateAndControls(T &x, Vectorcd &u, 
               const std::vector<VectorXd> &y) {

    std::cout << "[W] System::FlatToStateAndControls: unimplemented! Subclasses should override." << std::endl;
  }
}

#endif

