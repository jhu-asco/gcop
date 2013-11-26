#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "airmview.h"
#include "viewer.h"
#include "se3.h"
#include "utils.h"

using namespace gcop;
using namespace Eigen;

AirmView::AirmView(const Airm &sys, 
                   vector<MbsState> *xs) : 
  MbsView(sys, xs), hgv(sys.hrotor)
{
  this->geomViews.push_back(&hgv);
  for (int i = 0; i < 2; ++i) {
    //    SE3::Instance().rpyxyz2g(views[i].g0, Vector3d(0, 0, 0), Vector3d(0, 0, -sys.links[i+1].ds[2]/2));
    //    views[i].d = sys.links[i+1].ds;
    views[i].r = sys.links[i+1].ds[0]/2;
    views[i].h = sys.links[i+1].ds[2];

    this->geomViews.push_back(&views[i]);
  }
}
