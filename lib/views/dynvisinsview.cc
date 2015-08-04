#include "dynvisinsview.h"
#include <iostream>
#include "utils.h"

using namespace std;
using namespace gcop;

DynVisInsView::DynVisInsView(const DynVisIns &vi) : 
  View("DynVisual INS"), vi(vi), bodyView(body, (vector<Body3dState>*)&vi.xs)
{
  body.ds << .3, .1, .1;
}


DynVisInsView::~DynVisInsView()
{
}


void DynVisInsView::Render()
{
  ((SystemView<Body3dState,Vector6d>&)bodyView).Render();
  //  bodyView.Render();
  for (int j = 0; j < vi.ls.size(); ++j) {
    glPushMatrix(); 
    if (!vi.v) {
      glTranslated(vi.ls[j][0], vi.ls[j][1], vi.ls[j][2]);
    } else {
      int nx = (vi.optBias ? 15 : 9);
      int ind = nx*vi.xs.size() + 3*j;
      glTranslated(vi.v[ind], vi.v[ind + 1], vi.v[ind + 2]);
    }
    
    glutSolidSphere(.02, 5, 5);

    //if (vi.P.cols() >= 15 + (j+1)*3 && vi.P.rows() >= 15 + (j+1)*3)
    //       RenderEllipsoid(vi.P.block<3,3>(15 + j*3, 15 + j*3));
    
    glPopMatrix();    
  }
}
    
bool DynVisInsView::RenderFrame(int i)
{
  return bodyView.RenderFrame(i);
}


