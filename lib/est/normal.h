#ifndef GCOP_NORMAL_H
#define GCOP_NORMAL_H

#include <Eigen/Dense>
#include <vector>

namespace gcop {
  
  using namespace Eigen;
  using namespace std;

  class Normal {
  public:
    /**
     * n-dimensional normal distribution
     * @param n dimension
     */
    Normal(int n = 1);

    /**
     * n-dimensional normal distribution with mean mu and covariance P
     * @param mu mean
     * @param P covariance
     */    
    Normal(const VectorXd &mu, const MatrixXd &P);

    virtual ~Normal();

    /**
     * Compute likelihood of element x
     * @param x n-dimensional vector
     * @return likelihood of sample
     */
    double L(const VectorXd &x) const;


    /**
     * Sample from the distribution
     * @param x n-dimensional vector to be sampled
     * @return likelihood of sample
     */
    double Sample(VectorXd &x);    

    /**
     * Updates the Cholesky factor and the normalization constant
     * @return true if covariance is positive definite
     */
    bool Update();
    
    /**
     * Estimate the distribution using data xs and costs cs (optional)
     * @param xps data points and corresponding probabilities (should sum up to 1)
     */
    void Fit(const vector<pair<VectorXd, double> > xps);

    VectorXd mu;     ///< mean
    MatrixXd P;      ///< covariance
    
    double det;      ///< determinant
    MatrixXd Pinv;   ///< covariance inverse
    bool pd;         ///< covariance is positive-definite
    
    MatrixXd A;      ///< cholesky factor
    VectorXd rn;     ///< normal random vector

    double norm;     ///< normalizer

    int bd;          ///< force a block-diagonal structure with block dimension bd (0 by default means do not enforce)

    LLT<MatrixXd> llt; ///< LLT object to Cholesky
  };

}



#endif
