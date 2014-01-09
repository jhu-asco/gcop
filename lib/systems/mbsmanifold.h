#ifndef GCOP_MBSMANIFOLD_H
#define GCOP_MBSMANIFOLD_H

#include "manifold.h"
#include "mbsstate.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  // state dimension for nb-body system with 1-dof joints
#define MBS_DIM(nb) (12 + 2*(nb - 1))
#define MBS_DIM_FIXED(nb) (2*(nb - 1))
  
  class MbsManifold : public Manifold<MbsState> {
    
  public:    

    /**
     * @param nb number of bodies (including a base body)
     * @param fixed is the base body fixed (if fixed then it is treated as a "virtual" fixed body)
     */
    MbsManifold(int nb, bool fixed = false);
    
    void Lift(VectorXd &v,
              const MbsState &xa,
              const MbsState &xb);      

    void Retract(MbsState &xb, 
                 const MbsState &xa,
                 const VectorXd &v);

    void dtau(MatrixXd &M, const VectorXd &v);

    void Adtau(MatrixXd &M, const VectorXd &v);

    //    bool fixed;  ///< whether it has a fixed base

    bool cay;   ///< use Cayley instead of Exponential

  };
}


#endif
