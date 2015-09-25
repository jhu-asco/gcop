#include "dynvisinsview.h"
#include <iostream>
#include "utils.h"

using namespace std;
using namespace gcop;

DynVisInsView::DynVisInsView(const DynVisIns &vi) : 
  View("DynVisual INS"), vi(vi), bodyView(body)
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
  map<int, DynVisIns::Camera>::const_iterator camIter;
  for (camIter = vi.cams.begin(); camIter != vi.cams.end(); ++camIter) {
    const DynVisIns::Camera &cam = camIter->second;
    bodyView.Render(&cam.x);    
  }

  map<int, DynVisIns::Point>::const_iterator pntIter;
  int j = 0;
  for (pntIter = vi.pnts.begin(); pntIter != vi.pnts.end(); ++pntIter, ++j) {
    const DynVisIns::Point &p = pntIter->second;
    glPushMatrix(); 
    if (!vi.v) {
      glTranslated(p.l[0], p.l[1], p.l[2]);
    } else {
      int nx = (vi.optBias ? 15 : 9);
      int ind = nx*vi.cams.size() + 3*j;
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
  return false;
  //  return bodyView.RenderFrame(i);
}
