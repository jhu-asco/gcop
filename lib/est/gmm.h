#ifndef GCOP_GMM_H
#define GCOP_GMM_H

#include "normal.h"

namespace gcop {
  using namespace Eigen;
  using namespace std;
  
  /**
   * Gaussian mixture model
   *
   */
  class Gmm {
  public:
    /**
     * Construct a GMM with dimension n and k clusters
     * @param n dimension
     * @arapm k number of modes
     */
    Gmm(int n, int k = 1);

    virtual ~Gmm();

    bool Update();

    double L(const VectorXd &x) const;

    double Sample(VectorXd &x);
 
    void Init(const VectorXd &xlb, const VectorXd &xub);
    
    /**
     * Fit a GMM to data
     * @param xps data (pairs of vectors and weights), weights should add to one
     * @param a smoothing factor for updating the parameter v according to
     *         [ v_new = a*v_est + (1-a)*v_old ]; it set to 1 by default
     * @param iter maximum number of EM iterations (only applies for multiple mixtures)
     * @param S use additional noise matrix during EM for stability
     */
    void Fit(const vector<pair<VectorXd, double> > &xps, double a = 1,
             int iter = 50, const MatrixXd *S = 0);
    
    int k;              ///< number of mixture components
    vector<Normal> ns;  ///< normal distribution for each component   
    vector<double> ws;  ///< component weights
    vector<double> cdf; ///< CDF

    double tol; 
  };
}

#endif
