#ifndef GCOP_UUVVIEW_H
#define GCOP_UUVVIEW_H

#include "uuv.h"
#include "systemview.h"
#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#include "viewer.h"
#include "so3.h"


namespace gcop {

#include "utils.h"

  using namespace Eigen;
  
  template <int c = 6>
    class UuvView : public SystemView<UuvState, Matrix<double, c, 1> > {

    typedef Matrix<double, c, 1> Vectorcd;

  public:
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     */
    
    UuvView(const Uuv<c> &sys);
        
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    UuvView(const Uuv<c> &sys,
            vector<UuvState> *xs,
            vector<Vectorcd> *us = 0);
    
    virtual ~UuvView();
    
    
    virtual void Render(const UuvState *x,
                        const Vectorcd *u = 0);
    
    void Render(const vector<UuvState> *xs, 
                const vector<Vectorcd> *us = 0, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Uuv<c> &sys;
    GLUquadricObj *qobj;
    double dirSize;
    
    static void Transform(const Matrix4d &g);
  };
}


using namespace gcop;
using namespace Eigen;

template<int c>
UuvView<c>::UuvView(const Uuv<c> &sys) : 
  SystemView<UuvState, Matrix<double, c, 1> >("Uuv"), sys(sys)
{
  this->rgba[0] = .5;
  this->rgba[1] = .5;
  this->rgba[2] = .5;
  this->rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
  dirSize = -1;
}

template<int c>
UuvView<c>::UuvView(const Uuv<c> &sys,
                          vector<UuvState> *xs,
                          vector<Vectorcd> *us) : 
SystemView<UuvState, Matrix<double, c, 1> >("Uuv", xs, us), sys(sys)
{
  this->rgba[0] = 0.5;
  this->rgba[1] = 0.5;
  this->rgba[2] = 0.5;
  this->rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
}


template<int c>
UuvView<c>::~UuvView()
{
  free(qobj);
}

template<int c>
void UuvView<c>::Render(const UuvState *x,
                        const Vectorcd *u)
{
  //   glColor4f(1,0.5,0.5,0.5);
  
  glPushMatrix();
  Transform(x->g);
  glScaled(sys.ds[0], sys.ds[1], sys.ds[2]); 
  glutSolidCube(1);
  glPopMatrix();
}


template<int c>
void UuvView<c>::Render(const vector<UuvState> *xs, 
                        const vector<Vectorcd > *us, 
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
    ie = xs->size()-1;

  assert(is >= 0 && is <= xs->size()-1 && ie >= 0 && ie <= xs->size()-1);
  assert(is <= ie);

  glDisable(GL_LIGHTING);
  glLineWidth(this->lineWidth);
  glBegin(GL_LINE_STRIP);
  for (int i = is; i <= ie; i+=dit) {
    const UuvState &x = (*xs)[i];
    glVertex3d(x.g(0,3), x.g(1,3), x.g(2,3));
  }
  glEnd();
  glLineWidth(1);
  glEnable(GL_LIGHTING);
  
  if (rs) {
    for (int i = 0; i < xs->size(); i+=dis) {
      Render(&(*xs)[i]);
    }
  }

  if (dl)
    Render(&xs->back());
}

template<int c>
void UuvView<c>::Transform(const Matrix4d &g)
{
  const SO3 &so3 = SO3::Instance();
  glTranslated(g(0,3), g(1,3), g(2,3));
  
  Vector3d e;
  so3.log(e, g.topLeftCorner<3,3>());
  double n = e.norm();
  if (n > so3.tol) {
    e = e/n;
    glRotated(RAD2DEG(n), e[0], e[1], e[2]);
  }
}


#endif
