#ifndef GCOP_MBSMANIFOLD_H
#define GCOP_MBSMANIFOLD_H

#include "manifold.h"
#include "mbsstate.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  // state dimension for nb-body system
#define MBS_DIM(nb) (12 + 2*(nb - 1))
  
  class MbsManifold : public Manifold<MbsState> {
    
  public:    

    MbsManifold(int nb);
    
    void Lift(VectorXd &v,
              const MbsState &xa,
              const MbsState &xb);      

    void Retract(MbsState &xb, 
                 const MbsState &xa,
                 const VectorXd &v);

    void dtau(MatrixXd &M, const VectorXd &v);

    void Adtau(MatrixXd &M, const VectorXd &v);

  };
}


#endif
