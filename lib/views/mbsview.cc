#include "mbsview.h"

using namespace gcop;

MbsView::MbsView(const Mbs &sys,
                 vector<MbsState> *xs) : 
  SystemView<MbsState>("Mbs", xs), sys(sys)
{
  this->rgba[0] = 0.5;
  this->rgba[1] = 0.5;
  this->rgba[2] = 0.5;
  this->rgba[3] = 0;
  this->lineWidth = 2;
  this->qobj = gluNewQuadric();
}


MbsView::~MbsView()
{
  free(this->qobj);
}

void MbsView::Render(const MbsState &x)
{
  //   glColor4f(1,0.5,0.5,0.5);
  assert(geomViews.size() == sys.nb);
  
  for (int i = 0; i < sys.nb; ++i) {
    glPushMatrix();
    
    Geom3dView::Transform(x.gs[i]);
    
    // glScaled(sys.links[i].ds[0], sys.links[i].ds[1], sys.links[i].ds[2]); 
    assert(geomViews[i]);
    geomViews[i]->Render();//Geom();
    
    //    glutSolidCube(1);
    //   glutSolidSphere(.5, 10, 10);
    glPopMatrix();
  }
}


void MbsView::Render(const vector<MbsState> &xs, 
                     bool rs, 
                     int is, int ie,
                     int dis, int dit,
                     bool dl)
{
  Viewer::SetColor(this->rgba[0], this->rgba[1], this->rgba[2], this->rgba[3]);
  //  glColor4f(
  
  // set defaults
  if (is == -1)
    is = 0;
  if (ie == -1)
    ie = xs.size()-1;
  
  assert(is >= 0 && is <= xs.size()-1 && ie >= 0 && ie <= xs.size()-1);
  assert(is <= ie);

  glDisable(GL_LIGHTING);
  glLineWidth(this->lineWidth);
  glBegin(GL_LINE_STRIP);
  for (int i = is; i <= ie; i+=dit) {
    const MbsState &x = xs[i];
    glVertex3d(x.gs[0](0,3), x.gs[0](1,3), x.gs[0](2,3));
  }
  glEnd();
  glLineWidth(1);
  glEnable(GL_LIGHTING);
  
  if (rs) {
    for (int i = 0; i < xs.size(); i+=dis) {
      Render(xs[i]);
    }
  }

  if (dl)
    Render(xs.back());
}
