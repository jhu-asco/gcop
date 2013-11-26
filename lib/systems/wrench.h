#ifndef GCOP_WRENCH_H
#define GCOP_WRENCH_H

namespace gcop {
  
  using namespace Eigen;
  
  /**
   * An abstract wrench force function \f$ f:X\times U \rightarrow \mathbb{R}^6 \f$
   * depending on the state on control inputs, i.e. \f$ f(x,u) \f$ 
   *
   * This can be regarded as 6-dimensional wrench applied to the body
   * resulting from either control inputs, external forces, or both.
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <class T, int n, int c> class Wrench {
  public:

    typedef Matrix<double, 6, 1> Vector6d;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, c, c> Matrixcd;
    typedef Matrix<double, 6, c> Matrix6xcd;
    typedef Matrix<double, 6, n> Matrix6xnd;
            
    /**
     * Wrench force in 3d
     * @param dx whether it depends on the state x
     * @param du whether it depends on the controls u
     */
    Wrench(bool dx = true, bool du = true);

    
    /**
     * The force function 
     * @param f the resulting force
     * @param t time
     * @param x state
     * @param u controls
     * @param fx derivative w.r.t. state (optional)
     * @param fu derivative w.r.t. controls (optional)
     */
    virtual void F(Vector6d &f,
                   double t, const T &x, const Vectorcd &u,
                   Matrix6xnd *fx = 0, Matrix6xcd *fu = 0) = 0;
    
    bool dx;   ///< set to true if f depends on x
    bool du;   ///< set to true if f depends on u
  };
  
  template <class T, int n, int c> 
    Wrench<T, n, c>::Wrench(bool dx, bool du) :
    dx(dx), du(du) {
  } 
}

#endif
