#ifndef GCOP_RNLQSENSORCOST_H
#define GCOP_RNLQSENSORCOST_H

#include "lqsensorcost.h"
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
    int _ng = Dynamic,
    int _nz = Dynamic> class RnLqSensorCost :
    public LqSensorCost< Matrix<double, _nx, 1>, _nx, _nu, _np, _ng, Matrix<double, _nz, 1>, _nz> {
    
    typedef Matrix<double, _nx, 1> Vectornd;
    typedef Matrix<double, _nu, 1> Vectorcd;
    typedef Matrix<double, _np, 1> Vectormd;

    typedef Matrix<double, _nz, 1> Vectorrd;
    typedef Matrix<double, _nz, _nz> Matrixrd;
    typedef Matrix<double, _nz, _nx> Matrixrnd;
    typedef Matrix<double, _nz, _nu> Matrixrcd;
    typedef Matrix<double, _nz, _np> Matrixrmd;
    
  public:  
    /**
     * Linear-quadratic cost on a vector space. Use this constructor only for fixed-size 
     * initialization, i.e. RnLqCost<n,c>(tf, xf, ...)
     * @param diag are the R,S,P matrices diagonal? (true by default)
     */
    RnLqCost(System<Vectornd, _nx, _nu, _np> &sys, Manifold<Vectorrd, _nz> &Z, bool diag = true);
    
  };
  
  template <int _nx, int _nu, int _np, int _ng, int _nz> 
    RnLqSensorCost<_nx, _nu, _np, _ng, _nz>::RnLqCost(System<Vectornd, _nx, _nu, _np> &sys, 
                                                      Manifold<Vectorrd, _nz> &Z,
                                                      bool diag) : 
    LqSensorCost<Matrix<double, _nx, 1>, _nx, _nu, _np, _ng, Matrix<double, _nz, 1>, _nz>(sys, Z, diag) {    
    assert(_nx > 0);
    assert(_nu >= 0);
  }
}

#endif
