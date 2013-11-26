#ifndef GCOP_MODEL_H
#define GCOP_MODEL_H

#include <Eigen/Dense>
#include <cstring>

namespace gcop {

  /**
   * A general system model.
   * 
   * Author: Marin Kobilarov -- Copyright (C) 2005
   */
  class Model {
  public:    
    /**
     * System model that supports standard filtering on vector spaces, 
     * as well as constrained spaces that are for instance encountered
     * attitude filtering. In this case part of state represents a quaternion
     * which is constrained to have unit length, and so nc=1.
     *
     * @param nx state dimension
     * @param nz measurement dimension
     * @param nu input dimension (zero by default)
     * @Param nc number of constraints
     */
    Model(int nx, int nz, int nu = 0, int nc = 0);

    virtual ~Model();
    
    /**
     * Time update
     * @param xn next state
     * @param x current state
     * @return true if success
     */
    virtual bool f(Eigen::VectorXd &xn,
                   const Eigen::VectorXd &x);

    /**
     * Time update function
     * @param xn next state
     * @param x current state
     * @param u inputs
     * @return true if success
     */
    virtual bool f(Eigen::VectorXd &xn,
                   const Eigen::VectorXd &x,
                   const Eigen::VectorXd &u);
    
    /**
     * Measurement function
     * @param z measurement
     * @param x state
     * @return true if success
     */
    virtual bool h(Eigen::VectorXd &z, 
                   const Eigen::VectorXd &x);


    /**
     * Time-update function jacobian wrt state
     * @param A jacobian
     * @param x state
     * @return true if success
     */
    virtual bool dfdx(Eigen::MatrixXd &A,
                      const Eigen::VectorXd &x);    

    /**
     * Time-update function jacobian wrt state
     * @param A jacobian
     * @param x state
     * @param u inputs
     * @return true if success
     */
    virtual bool dfdx(Eigen::MatrixXd &A,
                      const Eigen::VectorXd &x,
                      const Eigen::VectorXd &u);

    /**
     * Time-update function jacobian wrt control
     * @param B jacobian
     * @param x state
     * @param u inputs
     * @return true if success
     */
    virtual bool dfdu(Eigen::MatrixXd &B,
                      const Eigen::VectorXd &x,
                      const Eigen::VectorXd &u);
        
    /**
     * Measurement function jacobian wrt state
     * @param H jacobian
     * @param x state
     * @return true if success
     */
    virtual bool dhdx(Eigen::MatrixXd &H, 
                      const Eigen::VectorXd &x);

    int nx;  ///< state dimension
    int nz;  ///< measurement dimension
    int nu;  ///< control input dimension
    int nc;  ///< number of constraints
    int nr;  ///< reduced dimension, i.e. nr=nx-nc
    
    bool quat;     ///< special flag denoting that this model has a state whose first 4 components correspond to a quaternion; note that in this case we must have that nr < nx, since there is at least one state constraint due to the quaternion belonging to \f$ S^3 \f$
    
    Eigen::MatrixXd Q;   ///< process noise
    Eigen::MatrixXd R;   ///< measurement noise
  };
}

#endif
