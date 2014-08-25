#include "ba.h"

using namespace gcop;
using namespace Eigen;
using namespace std;

Ba::Ba(Posegraph2d &pg) : 
  pg(pg), sys(), cost(sys, pg.ts.back(), pg), pddp(sys, cost, pg.ts, pg.gs, pg.us, pg.p)
{
  
}

