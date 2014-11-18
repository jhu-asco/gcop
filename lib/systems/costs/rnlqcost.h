#ifndef GCOP_RNLQCOST_H
#define GCOP_RNLQCOST_H

#include "lqcost.h"
#include "rn.h"
#include "system.h"

namespace gcop {
  
  using namespace Eigen;
  
  /**
   *  Linear-quadratic cost over a vector space \f$ \mathbb R^n \f$
   */
  template <int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic, 
    int _ng = Dynamic> class RnLqCost :
    public LqCost< Matrix<double, _nx, 1>, _nx, _nu, _np, _ng> {
    
    typedef Matrix<double, _nx, 1> Vectornd;
    typedef Matrix<double, _nu, 1> Vectorcd;
    typedef Matrix<double, _np, 1> Vectormd;

    
  public:  
    /**
     * Linear-quadratic cost on a vector space. Use this constructor only for fixed-size 
     * initialization, i.e. RnLqCost<n,c>(tf, xf, ...)
     * @param tf final time
     * @param xf desired final state
     * @param diag are the Q,Qf,R matrices diagonal? (true by default)
     */
    RnLqCost(System<Vectornd, _nx, _nu, _np> &sys, double tf, const Vectornd &xf, bool diag = true);
    
  };
  
  template <int _nx, int _nu, int _np, int _ng> 
    RnLqCost<_nx, _nu, _np, _ng>::RnLqCost(System<Vectornd, _nx, _nu, _np> &sys, 
                                           double tf, 
                                           const Matrix<double, _nx, 1> &xf, 
                                           bool diag) : 
    LqCost<Matrix<double, _nx, 1>, _nx, _nu, _np, _ng>(sys, tf, xf, diag) {    
    assert(_nx > 0);
    assert(_nu >= 0);
  }

}

#endif
