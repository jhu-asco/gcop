#ifndef GCOP_SYSTEM_H
#define GCOP_SYSTEM_H

#include <Eigen/Dense>
#include <vector>
#include <assert.h>
#include "manifold.h"
#include <iostream>

/*! \mainpage Robotic Planning and Control (GCOP) software library
 *
 * \section Intro
 *
 * This is a partial distribution containing only discrete optimal control functionality.
 * Several example systems are included such as a simple car model and a quadrotor model.
 *
 * In addition, a parameter-dependent discrete optimal control method is included for solving
 * sys id / adaptive control / parameter estimation problems. It is applied to a 
 * SLAM/bundle adjustment problem.
 * 
 * 
 * \section Installation
 * \subsection Build requirements
 *  g++, cmake, libeigen, opengl, glut
 *
 * \subsection Download
 * 
 *
 * \subsection Compilation
 * - To compile >: 
 * -- tar xfz gcop.tar.gz; cd gcop 
 * -- mkdir build; cd build; cmake ..; make
 * - To test an optimal control example>: cd bin; ./qrotortest, ./cartest, etc...
 * - To test a bundle adjustment example>: cd bin; ./batest
 *
 * \section API Reference
 * see docs/html
 *
 * \section Usage
 *
 *
 * \section Example
 *  see directory bin
 * \section Author
 *  Marin Kobilarov -- Copyright (C) 2005-2013
 * \section Keywords
 * discrete optimal control, discrete mechanics, motion planning, nonlinear control
 */

namespace gcop {
  
  using namespace Eigen;
  
  
  /**
   * Control system interface for dynamic systems on manifolds. 
   * Defines the discrete dynamics and any structural/inertial system parameters 
   * such as body dimensions, mass matrices, etc... 
   *
   * Subclasses should provide implementation of either a regular dynamics function
   * F or a parameter-dependent cost function Fp
   * (e.g. for sys id / adaptive control / optimal design problems)
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <typename T, typename Tu, int _n = Dynamic, int _c = Dynamic> class System {
    
  public:
  
  typedef Matrix<double, _n, 1> Vectornd;
  typedef Matrix<double, _c, 1> Vectorcd;
  typedef Matrix<double, _n, _n> Matrixnd;
  typedef Matrix<double, _n, _c> Matrixncd;
  typedef Matrix<double, _c, _n> Matrixcnd;
  typedef Matrix<double, _c, _c> Matrixcd;
  
  typedef Matrix<double, Dynamic, 1> Vectormd;
  typedef Matrix<double, Dynamic, Dynamic> Matrixmd;
  typedef Matrix<double, _n, Dynamic> Matrixnmd;
  typedef Matrix<double, Dynamic, _n> Matrixmnd;
  
  System(Manifold<T, _n> &X, Manifold<Tu, _c> &U);
  
  /**
   * Discrete dynamics update. The function computes a change in the state
   * v, regarded as a vector field which can then be used to update the new state 
   * x_b the old state x_a according to x_b = x_a + v (that is when the space is 
   * Euclidean). In general manifolds the update is x_b = Retract(x_a, v), i.e.
   * a retraction mapping is used to update x_b intrinsically.
   * @param v resulting change in state
   * @param t current time
   * @param xa current state
   * @param u current control
   * @param h current time-step
   * @param A jacobian of v(x)
   * @param B jacobian of v(u)
   */
  virtual double F(Vectornd &v, double t, const T& xa, 
                   const Tu &u, double h,
                   Matrixnd *A = 0, Matrixncd *B = 0);
  

  /**
   * Discrete dynamics update. The function computes a change in the state
   * v, regarded as a vector field which can then be used to update the new state 
   * x_b the old state x_a according to x_b = x_a + v (that is when the space is 
   * Euclidean). In general manifolds the update is x_b = Retract(x_a, v), i.e.
   * a retraction mapping is used to update x_b intrinsically.
   * @param v resulting change in state
   * @param t current time
   * @param xa current state
   * @param u current control
   * @param h current time-step
   * @param p static parameters
   * @param A jacobian of f(x)
   * @param B jacobian of f(u)
   * @param C jacobian of f(p)
   */
  virtual double F(Vectornd &v, double t, const T& xa,
                   const Tu &u, double h,
                   const Vectormd &p,
                   Matrixnd *A = 0, Matrixncd *B = 0, Matrixnmd *C = 0);
  
  virtual double Step(T& xb, double t, const T& xa,
                      const Tu &u, double h,
                      Matrixnd *A = 0, Matrixncd *B = 0);
  
  virtual double Step(T& xb, double t, const T& xa,
                      const Tu &u, double h,
                      const Vectormd &p,
                      Matrixnd *A = 0, Matrixncd *B = 0, Matrixnmd *C = 0);
  
  
  /**
   * Reconstruct the full state. Some states contain redundant parameters 
   * for efficiency and it is useful to reconstruct them to maintain a 
   * consisten state.
   * @param x state
   * @param h time-step
   */
  virtual void Rec(T &x, double h) {};
  
  Manifold<T, _n> &X;   ///< state manifold 
  Manifold<Tu, _c> &U;   ///< control manifold
  
  int n;   ///< state dimension
  int c;   ///< control dimension
  
  };
  
  
  
  template <typename T, typename Tu, int _n, int _c> 
    System<T, Tu, _n, _c>::System(Manifold<T, _n> &X, Manifold<Tu, _c> &U) : 
    X(X), U(U), n(X.n), c(U.n) {
    assert(n > 0);
    assert(c >= 0);
  }

  
  template <typename T, typename Tu, int _n, int _c> 
    double System<T, Tu, _n, _c>::F(Vectornd &v, double t, const T &xa, 
                                    const Tu &u, double h,
                                    Matrix<double, _n, _n> *A, Matrix<double, _n, _c> *B) {
    std::cout << "[W] System::F: unimplemented!" << std::endl;
    return 0;
  }
  
  template <typename T, typename Tu, int _n, int _c> 
    double System<T, Tu, _n, _c>::F(Vectornd &v, double t, const T& xa,
                                     const Tu &u, double h,
                                     const VectorXd &p,
                                     Matrix<double, _n, _n> *A, Matrix<double, _n, _c> *B, 
                                     Matrix<double, _n, Dynamic> *C) {
    std::cout << "[W] System::F: unimplemented!" << std::endl;
    return 0;
  }

  template <typename T, typename Tu, int _n, int _c> 
    double System<T, Tu, _n, _c>::Step(T &xb, double t, const T &xa,
                                       const Tu &u, double h,
                                       Matrix<double, _n, _n> *A, Matrix<double, _n, _c> *B) {
    Vectornd v;
    if (_n == Dynamic)
      v.resize(n);

    double d = F(v, t, xa, u, h, A, B);
    X.Retract(xb, xa, v); // in R^n this is just xb = xa + v
    Rec(xb, h);
    return d;
  }
  
  template <typename T, typename Tu, int _n, int _c> 
    double System<T, Tu, _n, _c>::Step(T& xb, double t, const T& xa,
                                       const Tu &u, double h,
                                       const VectorXd &p,
                                       Matrix<double, _n, _n> *A, Matrix<double, _n, _c> *B, 
                                       Matrix<double, _n, Dynamic> *C) {
    Vectornd v;
    if (_n == Dynamic)
      v.resize(n);

    double d = F(v, t, xa, u, h, p, A, B, C);
    X.Retract(xb, xa, v); // in R^n this is just xb = xa + v
    Rec(xb, h);
    return d;
  }
}

#endif
