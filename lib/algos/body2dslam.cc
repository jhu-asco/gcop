#include "body2dslam.h"

using namespace gcop;
using namespace Eigen;
using namespace std;

Body2dSlam::Body2dSlam(Body2dGraph &pg) : 
  pg(pg), cost(pg.ts.back(), pg)
{
  pdmoc = new PDmoc<M3V3d, 6, 3>(pg.sys, cost, pg.ts, pg.xs, pg.us, pg.p, 2*pg.extforce);
}
