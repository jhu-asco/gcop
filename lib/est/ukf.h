#ifndef GCOP_UKF_H
#define GCOP_UKF_H

#include "filter.h"

namespace gcop {
    
  /**
   * General Unscented Kalman Filter. This is the additive version, 
   * i.e. the state is not agumented with noise terms.
   * 
   * Author: Marin Kobilarov -- Copyright (C) 2005
   */
  class UKF : public Filter {
  public:
    /**
     * Basic Kalman Filter
     * @param model system model
     */
    UKF(Model &model);

    virtual ~UKF();
        
    /**
     * Prediction step     
     * @param u control (optional - no control by default)
     * @param cov whether to update the covariance as well (true by default)
     * @return true on success
     */
    virtual bool Predict(const Eigen::VectorXd *u = 0,
                         bool cov = true);
    
    /**
     * Update step
     * @param z measurement
     * @param cov whether to update the covariance as well (true by default)
     * @return true on success
     */
    virtual bool Update(const Eigen::VectorXd &z, 
                        bool cov = true);

  protected:
    
    /**
     * Compute 2*L+1 sigma points given mean x and spread A
     * @param Xs pointer to an array of 2*L+1 vectors
     * @param x mean state
     * @param A the square root of the covariance 
     */
    void Points(Eigen::VectorXd **Xs,
                const Eigen::VectorXd &x,
                const Eigen::MatrixXd &A);
    
    double a;      ///< UKF \f$\alpha\f$-param (a=1e-3 by default)
    double k;      ///< UKF \f$\kappa\f$-param (k=0 by default)
    double b;      ///< UKF \f$\beta\f$-param (b=2 by default)

  protected:
    int L;         ///< number of sigma points (this is set internally to 2*model.nr + 1
    double l;      ///< the variable \f$ l = a^2(L+k)-L \f$

    Eigen::VectorXd** Xs;   ///< state sigma points
    Eigen::VectorXd** Xps;  ///< predicted state sigma points
    Eigen::VectorXd** Zs;   ///< measurement sigma points
    
    Eigen::VectorXd Ws;  ///< points weights
    Eigen::VectorXd Wc;  ///< cov update weights

    Eigen::MatrixXd Pzz; ///< internally used covariances
    Eigen::MatrixXd Pxz; ///< internally used covariances

    Eigen::MatrixXd A;   ///< cholesky factor
  };
}

#endif
