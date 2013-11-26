#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "chainview.h"
#include "viewer.h"
#include "so3.h"
#include "utils.h"
#include "boxview.h"

using namespace gcop;
using namespace Eigen;

ChainView::ChainView(const Chain &sys,
                     vector<MbsState> *xs) : 
  MbsView(sys, xs)
{
  for (int i = 0; i < 3; ++i) {
    views[i].d = sys.links[i].ds;
    this->geomViews.push_back(&views[i]);
  }
}
