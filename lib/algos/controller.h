#ifndef GCOP_CONTROLLER_H
#define GCOP_CONTROLLER_H

#include <Eigen/Dense>

namespace gcop {
  
  using namespace Eigen;
  
  /**
   * Base controller. Its simple purpose is to take a state x
   * and time t and produces the control u.
   *
   * Author: Marin Kobilarov -- Copyright (C) 2006-2013
   */
  template <typename Tx, typename Tu, typename Ts = VectorXd, typename Tc = VectorXd> class Controller {
  public:
  /**
   * Feedback controller of the form u = k(t, x)
   * @param u controls
   * @param t time
   * @param x state
   * @return true if arguments are valid, false if e.g. x is infeasible
   */
  virtual bool Set(Tu &u, double t, const Tx &x) = 0;
  
  /**
   * Set controller parameters
   * @return true if arguments are valid, false if e.g. s is infeasible
   */
  virtual bool SetParams(const Ts &s) {
    this->s = s;
    return true;
  }
  
  /**
   * Set context
   * @return true if arguments are valid, false if e.g. s is infeasible
   */
  virtual bool SetContext(const Tc &c) {
    this->c = c;
    return true;
  }
  
  Ts s; ///< controller parameters
  Tc c; ///< controller context
  };  
};

#endif
