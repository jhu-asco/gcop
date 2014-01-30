#ifndef GCOP_FORCE_H
#define GCOP_FORCE_H

#include <Eigen/Dense>
#include <vector>
#include <assert.h>
#include <iostream>

namespace gcop {
  
  using namespace Eigen;
  
  
  /**
   * Control force interface for dynamic forces on manifolds. 
   * Defines a force f (in R^nf) as a function of 
   * state x (of an n-dimensional manifold) and force u (of R^c)
   *
   * Subclasses should provide implementation of the function Setp
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <typename Tx, int _l = Dynamic, int _n = Dynamic, int _c = Dynamic> class Force {
    
  public:
  
  typedef Matrix<double, _l, 1> Vectorld;
  typedef Matrix<double, _n, 1> Vectornd;
  typedef Matrix<double, _c, 1> Vectorcd;
  typedef Matrix<double, _l, _n> Matrixlnd;
  typedef Matrix<double, _l, _c> Matrixlcd;
  
  typedef Matrix<double, _l, Dynamic> MatrixlXd;
  
  virtual void Set(Vectorld &f, double t, const Tx &x, 
                   const Vectorcd &u, double h, const VectorXd *p = 0,
                   Matrixlnd *A = 0, Matrixlcd *B = 0, MatrixlXd *C = 0);
  };
  
  template <typename Tx, int _l, int _n, int _c> 
    void Force<Tx, _l, _n, _c>::Set(Vectorld &f, double t, const Tx &x, 
                                    const Vectorcd &u, double h, const VectorXd *p,
                                    Matrix<double, _l, _n> *A, Matrix<double, _l, _c> *B,
                                    MatrixlXd *C) {
    std::cout << "[W] Force::Set: unimplemented!" << std::endl;
  }  
}

#endif
