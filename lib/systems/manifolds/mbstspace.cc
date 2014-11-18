#include "mbstspace.h"

using namespace gcop;

MbsTspace::MbsTspace(int nb) : Manifold<MbsState>(nb + 5) {
}


void MbsTspace::Lift(VectorXd &v, 
                     const MbsState &xa,
                     const MbsState &xb) {
  
  v.head(6) = xb.vs[0] - xa.vs[0];
  v.tail(xa.dr.size()) = xb.dr - xa.dr;
}


void MbsTspace::Retract(MbsState &xb, 
                        const MbsState &xa,
                        const VectorXd &v) {
  xb.vs[0] = xa.vs[0] + v.head(6);
  xb.dr = xa.dr + v.tail(xa.dr.size());
}
