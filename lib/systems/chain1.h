#ifndef GCOP_CHAIN1_H
#define GCOP_CHAIN1_H

#include "mbs.h"

namespace gcop {
  
  class Chain1 : public Mbs {
  public:
    Chain1();    

    void Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
               MatrixXd *A = 0, MatrixXd *B = 0);
  };
}

#endif
