#ifndef GCOP_KINBODY3DTRACKVIEW_H
#define GCOP_KINBODY3DTRACKVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include "kinbody3dtrack.h"
#include "view.h"


namespace gcop {

  /**
   * Rigid body view. Supports rendering either a single state
   * or a whole trajectory of states.
   */
  template <int _nu = 6>
  class Kinbody3dTrackView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    Kinbody3dTrackView(const Kinbody3dTrack<_nu> &pg);

    virtual ~Kinbody3dTrackView();

    virtual void Render();
    
    virtual bool RenderFrame(int i);
  
    const Kinbody3dTrack<_nu> &pg;

    float rgba[4];

    bool drawLandmarks;
    bool drawForces;
    double forceScale;


    GLUquadricObj *qobj;   
  };

template <int _nu>
Kinbody3dTrackView<_nu>::Kinbody3dTrackView(const Kinbody3dTrack<_nu> &pg) : 
  View("Kinbody3dtrack View"), pg(pg)
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

template <int _nu>
Kinbody3dTrackView<_nu>::~Kinbody3dTrackView()
{
  free(qobj);
}


template <int _nu>
void Kinbody3dTrackView<_nu>::Render()
{
  RenderFrame(0);
}


template <int _nu>
bool Kinbody3dTrackView<_nu>::RenderFrame(int i)
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
      f.head<3>() = forceScale*pg.xs[k].block(0,0,3,3)*(pg.sys.Bu*pg.us[k]).tail(3);
      glPushMatrix();
      glTranslated(pg.xs[k](0,3), pg.xs[k](1,3), pg.xs[k](2,3)); 
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
        const Vector3d &x = pg.xs[k].block(0,3,3,1);
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

}

#endif
