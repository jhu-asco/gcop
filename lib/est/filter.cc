#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include "filter.h"

using namespace gcop;
using namespace Eigen;


Filter::Filter(Model &model) : 
  model(model),
  x(model.nx),
  P(model.nr, model.nr)
{
  Init();
}

void Filter::Init()
{
  chi = 0;
  chiGate = -1;
}


Filter::~Filter()
{
}
