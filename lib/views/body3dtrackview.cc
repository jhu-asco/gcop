#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "viewer.h"
#include "utils.h"
#include "body3dtrackview.h"

using namespace gcop;

Body3dTrackView::Body3dTrackView(const Body3dTrack &pg) : 
  View("Body3dtrack View"), pg(pg)
{
  rgba[0] = 1;
  rgba[0] = 0;
  rgba[0] = 0;
  rgba[0] = .5;  
  drawLandmarks = true;
  drawForces = true;
  forceScale = 10;
  qobj = gluNewQuadric();
}


Body3dTrackView::~Body3dTrackView()
{
  free(qobj);
}


void Body3dTrackView::Render()
{
  RenderFrame(0);
}


bool Body3dTrackView::RenderFrame(int i)
{  

  Viewer::SetColor(.5, .5, .5, 1);  

  gluCylinder(qobj, pg.r + pg.w/2, pg.r + pg.w/2, .1, 20, 10);

  gluCylinder(qobj, pg.r - pg.w/2, pg.r - pg.w/2, .1, 20, 10);

  // red: noisy/estimated
  Viewer::SetColor(rgba[0], rgba[1], rgba[2], rgba[3]);

  if (drawForces) {

    glDisable(GL_LIGHTING);
    Vector3d f;
    for (int k =0; k < pg.us.size(); ++k) {
      f.head<3>() = forceScale*pg.xs[k].R*(pg.us[k]).tail<3>();
      glPushMatrix();
      glTranslated(pg.xs[k].p(0), pg.xs[k].p(1), pg.xs[k].p(2)); 
      Viewer::DrawArrow(f.data(),qobj);
      glPopMatrix();
    }
    glEnable(GL_LIGHTING);
  }

  if (drawLandmarks) {

    
    for (int l = 0; l < pg.ls.size(); ++l) {    

      if (pg.observed[l])
        Viewer::SetColor(0, 1, 0, 1);  
      else
        Viewer::SetColor(.5, .5, .5, 1);  

      glPushMatrix();
      glTranslated( pg.ls[l][0], pg.ls[l][1], pg.ls[l][2]);

      gluCylinder(qobj, pg.pr, pg.pr, 1, 20, 10);

      //      glutSolidSphere(pg.pr, 10, 10);
      
      /*
      glPushMatrix();
      glTranslated(.05, 0, 0);
      char str[10];
      sprintf(str, "%d", l);
      Viewer::DrawText(str, .05);
      glPopMatrix();
      */
      glPopMatrix();
    }

    // draw visible
    int nvf = (pg.p.size() - 3*pg.extforce)/3;
    int i0 = 3*pg.extforce;
        
    Viewer::SetColor(0, 0, 1, 1);  
    //    cout << "ALL" << endl;
    for (int l = 0; l < nvf; ++l) {    
      //      cout << pg.p[i0 + 2*l] << " " << pg.p[i0 + 2*l + 1] << endl;
      glPushMatrix();
      glTranslated( pg.p[i0 + 3*l], pg.p[i0 + 3*l + 1], pg.p[i0 + 3*l + 2]);
      glutSolidSphere(pg.pr, 10, 10);
    
      /*
      glPushMatrix();
      glTranslated(.05, 0, 0);
      char str[10];
      sprintf(str, "%d", l);
      Viewer::DrawText(str, .05);
      glPopMatrix();
      */
      glPopMatrix();
    }   
    
    glDisable(GL_LIGHTING);
    glLineWidth(1);
    glBegin(GL_LINES);
    //    Viewer::SetColor(0, 0, 1, 0);  
    glColor3d(.5,.5,1);
    for (int l = 0; l < nvf; ++l) {
      const vector< pair<int,Vector3d> > &J = pg.Js[pg.pis[l]];
      for (int j = 0; j < J.size(); ++j) {
        int k = J[j].first;
        const Vector3d &x = pg.xs[k].p;
        glVertex3d(pg.p(3*l + i0), pg.p(3*l + i0 + 1), pg.p(3*l + i0 + 2));
        glVertex3d(x[0], x[1], x[2]);
      }
      }
    glEnd();
    glLineWidth(1);
    glEnable(GL_LIGHTING);

  }

  return false;
}

