#ifndef GCOP_CONTROLLER_H
#define GCOP_CONTROLLER_H

namespace gcop {

  /**
   * Base controller. Its simple purpose is to take a state x
   * and time t and produces the control u.
   *
   * Author: Marin Kobilarov -- Copyright (C) 2006-2013
   */
  template <typename Tx, typename Tu> class Controller {
  public:
    /**
     * Feedback controller of the form u = k(t, x)
     * @param u controls
     * @param t time
     * @param x state
     */
    virtual void Set(Tu &u, double t, const Tx &x) = 0;
    
  };  
};

#endif
