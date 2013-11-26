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
    
    void Fit(const vector<pair<VectorXd, double> > &xps,
             int iter = 50, const MatrixXd *S = 0);
    
    int k;              ///< number of mixture components
    vector<Normal> ns;  ///< normal distribution for each component   
    vector<double> ws;  ///< component weights
    vector<double> cdf; ///< CDF

    double tol; 
  };
}

#endif
