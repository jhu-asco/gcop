// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_GP_H
#define GCOP_GP_H

#include <Eigen/Dense>

namespace gcop {
  class GP {
  public:
    /**
     * Initialize a GP with dimension d and number of points n
     * @param d dimension
     * @param n number of points (optional)
     */
    GP(int d, int n = 0);
    
    virtual ~GP();
    
    /**
     * Generate random points
     */ 
    virtual void Sample();
    
    /**
     * Add a new data point
     * @param x data vector
     * @param f value
     * @param true if OK
     */
    bool Add(const Eigen::VectorXd &x,
             double f);    
    
    /**
     * Train using current data
     */
    void Train();
    
    /**
     * Train GP using a given dataset (xs, fs)
     * @param d-n matrix of data vectors
     * @param n-vector of values
     */
    void Train(const Eigen::MatrixXd &Xs, 
               const Eigen::VectorXd &fs);
    
    /**
     * Predict value at point x
     * @param x point
     * @param s pointer to predicted covariance (optional)
     * @return predicted mean
     */
    double Predict(const Eigen::VectorXd &x, 
                   double *s = 0) const;
    
    /**
     * Square exponential kernel
     * @param xa first point
     * @param xb second point
     * @return correlation
     */
    double SqExp(const Eigen::VectorXd &xa, 
                 const Eigen::VectorXd &xb) const;
    
    /**
     * Loglikelihood
     * @param dll derivative of log-liklihood w.r. to l and s
     * @return log-likelihood
     */
    double LogL(double dll[2] = 0);

    /**
     * Probability of improvement over a given value fmin
     * @aram x data vector
     * @param fmin baseline
     * @return probability of improvement
     */
    double PI(const Eigen::VectorXd &x, double fmin) const;

    /**
     * Optimize GP parameters. This is currently done naively
     * using a grid enumeration over l and s
     */
    double OptParams();
       
    int d;  ///< dimension
    int n;  ///< number of data points
    
    Eigen::MatrixXd Xs;   ///< data points
    Eigen::VectorXd fs;   ///< values
    
    Eigen::MatrixXd K;  
    Eigen::MatrixXd Ki;
    
    Eigen::MatrixXd L;
    
    Eigen::VectorXd a;
    
    double l;
    double s;
    
    double sigma;
    
    bool cf;   ///< propagate cholesky factor L rather than K^{-1}

    bool eps;  ///< prohibit adding points that are eps-close in L_2 to existing data 
    
  };
}

#endif
