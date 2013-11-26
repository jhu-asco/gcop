#ifndef GCOP_RNLQCOST_H
#define GCOP_RNLQCOST_H

#include "lqcost.h"
#include "rn.h"

namespace gcop {
  
  using namespace Eigen;
  
  /**
   *  Linear-quadratic cost over a vector space \f$ \mathbb R^n \f$
   */
  template <int _n = Dynamic, int _c = Dynamic> class RnLqCost :
    public LqCost< Matrix<double, _n, 1>, _n, _c> {
    
    typedef Matrix<double, _n, 1> Vectornd;
    typedef Matrix<double, _c, 1> Vectorcd;

    
  public:  
    /**
     * Linear-quadratic cost on a vector space. Use this constructor only for fixed-size 
     * initialization, i.e. RnLqCost<n,c>(tf, xf, ...)
     * @param tf final time
     * @param xf desired final state
     * @param diag are the Q,Qf,R matrices diagonal? (true by default)
     */
    RnLqCost(double tf, const Vectornd &xf, bool diag = true);

    /**
     * Linear-quadratic cost on a vector space. Use this constructor for dynamic-size 
     * initialization, i.e. RnLqCost<>(n, c, tf, xf, ...)
     * @param X state space
     * @param U control space
     * @param tf final time
     * @param xf desired final state
     * @param diag are the Q,Qf,R matrices diagonal? (true by default)
     */
    RnLqCost(Rn<_n> &X, Rn<_c> &U, double tf, const Vectornd &xf, bool diag = true);
    
  };
  
  template <int _n, int _c> 
    RnLqCost<_n, _c>::RnLqCost(double tf, 
                               const Matrix<double, _n, 1> &xf, 
                               bool diag) : 
    LqCost<Matrix<double, _n, 1>, _n, _c>(Rn<_n>::Instance(), Rn<_c>::Instance(), tf, xf, diag) {    
    assert(_n > 0);
    assert(_c >= 0);
  }
  
  template <int _n, int _c> 
    RnLqCost<_n, _c>::RnLqCost(Rn<_n> &X, 
                               Rn<_c> &U, 
                               double tf, 
                               const Matrix<double, _n, 1> &xf, 
                               bool diag) : 
    LqCost<Matrix<double, _n, 1>, _n, _c>(X, U, tf, xf, diag) {    
    assert(X.n > 0);
    assert(U.c >= 0);
  }

}

#endif
