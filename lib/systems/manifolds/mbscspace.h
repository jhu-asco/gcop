#ifndef GCOP_MBSCSPACE_H
#define GCOP_MBSCSPACE_H

#include "manifold.h"
#include "mbsstate.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  // state dimension for nb-body system

  class MbsCspace : public Manifold<MbsState> {
  public:

    MbsCspace(int nb);

    void Lift(Vectornd &v,
              const MbsState &xa,
              const MbsState &xb);      

    void Retract(MbsState &xb, 
                 const MbsState &xa,
                 const VectorXd &v);
  };
}


#endif
