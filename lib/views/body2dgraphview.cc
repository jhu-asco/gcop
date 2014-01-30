#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "viewer.h"
#include "utils.h"
#include "body2dgraphview.h"

using namespace gcop;

Body2dGraphView::Body2dGraphView(const Body2dGraph &pg) : 
  View("Body2dgraph View"), pg(pg)
{
  rgba[0] = 1;
  rgba[0] = 0;
  rgba[0] = 0;
  rgba[0] = .5;  
  drawLandmarks = true;
}


void Body2dGraphView::Render()
{
  RenderFrame(0);
}


bool Body2dGraphView::RenderFrame(int i)
{  
  // red: noisy/estimated

  glDisable(GL_LIGHTING);
  glColor4d(rgba[0], rgba[1], rgba[2], rgba[3]);
  //  Viewer::SetColor(rgba[0], rgba[1], rgba[2], rgba[3]);  
  Vector3d f;
  for (int k =0; k < pg.us.size(); ++k) {
    f[2] = 0;
    f.head<2>() = pg.xs[k].first.topLeftCorner<2,2>()*pg.us[k].tail<2>();
    glPushMatrix();
    glTranslated(pg.xs[k].first(0,2), pg.xs[k].first(1,2), 0); 
    Viewer::DrawArrow(f.data());
    glPopMatrix();
  }
  glEnable(GL_LIGHTING);


  if (drawLandmarks) {
    int nf = (pg.p.size() - 2*pg.extforce)/2;
    int i0 = 2*pg.extforce;
    
    for (int l = 0; l < nf; ++l) {    
      glPushMatrix();
      glTranslated(pg.p[i0 + 2*l], pg.p[i0 + 2*l+1], 0);
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
    for (int l = 0; l < nf; ++l) {
      const vector< pair<int,Vector2d> > &J = pg.Js[l];
      assert(J.size());
      for (int j = 0; j < J.size(); ++j) {
        int k = J[j].first;
        const Vector2d &x = pg.xs[k].first.block<2,1>(0,2);
        glVertex3d(pg.p(2*l + i0), pg.p(2*l + i0 + 1), 0);
        glVertex3d(x[0], x[1], 0);
      }
    }
    glEnd();
    glLineWidth(1);
    glEnable(GL_LIGHTING);
  }

  return false;
}

