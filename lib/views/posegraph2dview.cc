#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "viewer.h"
#include "utils.h"
#include "posegraph2dview.h"

using namespace gcop;

Posegraph2dView::Posegraph2dView(const Posegraph2d &pg) : 
  View("Posegraph2d View"), pg(pg)
{
  rgba[0] = 1;
  rgba[0] = 0;
  rgba[0] = 0;
  rgba[0] = .5;  
}


void Posegraph2dView::Render()
{
  RenderFrame(0);
}


bool Posegraph2dView::RenderFrame(int i)
{
  
  // red: noisy/estimated
  Viewer::SetColor(rgba[0], rgba[1], rgba[2], rgba[3]);  

  for (int l = 0; l < pg.p.size()/2; ++l) {    
    glPushMatrix();
    glTranslated(pg.p[2*l], pg.p[2*l+1], 0);
    glutSolidSphere(.05, 5, 5);

    glPushMatrix();
    glTranslated(.05, 0, 0);
    char str[10];
    sprintf(str, "%d", l);
    Viewer::DrawText(str, .05);
    glPopMatrix();

    glPopMatrix();
  }

  glDisable(GL_LIGHTING);
  glLineWidth(1);
  glBegin(GL_LINES);
  int nf = pg.p.size()/2;
  for (int l = 0; l < nf; ++l) {
    const vector< pair<int,Vector2d> > &J = pg.Js[l];
    assert(J.size());
    for (int j = 0; j < J.size(); ++j) {
      int k = J[j].first;
      const Vector2d &x = pg.gs[k].block<2,1>(0,2);
      glVertex3d(pg.p(2*l), pg.p(2*l +1), 0);
      glVertex3d(x[0], x[1], 0);
    }
  }
  glEnd();
  glLineWidth(1);
  glEnable(GL_LIGHTING);

  return false;
}

