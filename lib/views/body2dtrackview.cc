#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "viewer.h"
#include "utils.h"
#include "body2dtrackview.h"

using namespace gcop;

Body2dTrackView::Body2dTrackView(const Body2dTrack &pg) : 
  View("Body2dtrack View"), pg(pg)
{
  rgba[0] = 1;
  rgba[0] = 0;
  rgba[0] = 0;
  rgba[0] = .5;  
  drawLandmarks = true;
  drawForces = true;
  forceScale = 1;
  qobj = gluNewQuadric();
}


Body2dTrackView::~Body2dTrackView()
{
  free(qobj);
}


void Body2dTrackView::Render()
{
  RenderFrame(0);
}


bool Body2dTrackView::RenderFrame(int i)
{  

  Viewer::SetColor(.5, .5, .5, 1);  

  gluCylinder(qobj, pg.r + pg.w/2, pg.r + pg.w/2, .1, 20, 10);

  gluCylinder(qobj, pg.r - pg.w/2, pg.r - pg.w/2, .1, 20, 10);

  // red: noisy/estimated
  //  Viewer::SetColor(rgba[0], rgba[1], rgba[2], rgba[3]);
  glColor4fv(rgba);

  if (drawForces) {

    glDisable(GL_LIGHTING);
    Vector3d f;
    for (int k =0; k < pg.us.size(); ++k) {
      f[2] = 0;
      f.head<2>() = forceScale*pg.xs[k].first.topLeftCorner<2,2>()*(pg.us[k]-pg.uos[k]).tail<2>();
      glPushMatrix();
      glTranslated(pg.xs[k].first(0,2), pg.xs[k].first(1,2), 0); 
      Viewer::DrawArrow(f.data(),qobj);
      glPopMatrix();
    }
    glEnable(GL_LIGHTING);
  }

  if (drawLandmarks) {
    for (int l = 0; l < pg.ls.size(); ++l) {    

      if (pg.observed[l]) {
        //        Viewer::SetColor(0, 1, 0, 1);  
        glColor3f(0,1,0);
      } else {
        //        Viewer::SetColor(.5, .5, .5, 1);  
        glColor3f(.5,.5,.5);
      }

      glPushMatrix();
      glTranslated( pg.ls[l][0], pg.ls[l][1], 0);

      //      gluCylinder(qobj, pg.pr, pg.pr, 1, 20, 10);

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

    // draw visible
    int nvf = (pg.p.size() - 2*pg.extforce)/2;
    int i0 = 2*pg.extforce;
        
    //    Viewer::SetColor(0, 0, 1, 1);  
    glColor3f(0,0,1);

    //    cout << "ALL" << endl;
    for (int l = 0; l < nvf; ++l) {    
      //      cout << pg.p[i0 + 2*l] << " " << pg.p[i0 + 2*l + 1] << endl;
      glPushMatrix();
      glTranslated( pg.p[i0 + 2*l], pg.p[i0 + 2*l + 1], 0);
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
      const vector< pair<int,Vector2d> > &J = pg.Js[pg.pis[l]];
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

