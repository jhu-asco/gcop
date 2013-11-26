#ifndef GCOP_CE_H
#define GCOP_CE_H

#include <vector>
#include "gmm.h"

namespace gcop {

  using namespace Eigen;
  using namespace std;
  
  /**
   * Cross-entropy optimization method. One iteration of the algorithm typically includes the following steps:
   *
   * 0. Initialize the GMM distribution, select number of samples N, and quantile rho
   *
   * 1. for 1:N, z = Sample(), c = Cost(z), AddSample(z,c), end
   *
   * 2. Select() selects the top rho-quantile by sorting the samples using their costs
   *
   * 3. Fit() estimates the distribution
   *
   * 4. Best() gives the best sample or just take the first sample in the list of samples as the solution
   *
   * 5. Reset() and goto 1.
   *
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Ce {
  public:    
    
    /**
     * Initialize n-dimensional CE with k modes
     * @param n parameter space dimension
     * @param k number of mixture components
     * @param S n-n matrix of random noise to be added to the covariance after GMM estimation; acts both as regularizing factor and also for improving local exploration. This should typically be set to a small value, e.g. a diagonal matrix with entries = 0.01
     */
    Ce(int n,
       int k = 1,
       const MatrixXd *S = 0);
    
    virtual ~Ce();

    /**
     * Clear samples
     */
    void Reset();

    /**
     * Add a sample to list of samples
     * @param z sample
     * @param c cost
     */
    virtual void AddSample(const VectorXd &z, double c);

    /**
     * Select the top quantile. Resizes samples to size N*rho containing the best samples
     */
    void Select();
    
    /**
     * Fit GMM to data
     */
    bool Fit();
    
    /**
     * Draw a sample from the GMM
     * @param z vector to be sampled
     * @return likelihood
     */
    virtual double Sample(VectorXd &z);

    /**
     * Get the best sample. Equivalent to accessing zps[0].first. This assumes
     * Select() was called to sort the samples.
     * @return the sample with lowest cost
     */
    const VectorXd& Best();

    int n;          ///< parameter space dimension
      
    Gmm gmm;        ///< underlying Gaussian-mixture-model

    MatrixXd S;   ///< extra noise

    vector<pair<VectorXd, double> > zps;  ///< samples

    vector<double> cs;    ///< costs

    double rho;     ///< quantile

    double alpha;   ///< update factor

    bool mras;      ///< whether to use MRAS version of algorithm

    double b;       ///< MRAS rate
        
    //    Sample bsample;  ///< best sample
        
  };
}


#endif
