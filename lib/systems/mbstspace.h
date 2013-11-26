#ifndef GCOP_MBSTSPACE_H
#define GCOP_MBSTSPACE_H

#include "manifold.h"
#include "mbsstate.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  // state dimension for nb-body system

  class MbsTspace : public Manifold<MbsState> {
    
  public:

    MbsTspace(int nb);
    
    void Lift(VectorXd &v,
              const MbsState &xa,
              const MbsState &xb);      

    void Retract(MbsState &xb, 
                 const MbsState &xa,
                 const VectorXd &v);
  };
}


#endif
