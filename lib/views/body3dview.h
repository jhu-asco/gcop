#ifndef GCOP_BODY3DVIEW_H
#define GCOP_BODY3DVIEW_H

#include "body3d.h"
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
    class Body3dView : public SystemView<Body3dState, Matrix<double, c, 1> > {

    typedef Matrix<double, c, 1> Vectorcd;

  public:
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     */
    
    Body3dView(const Body3d<c> &sys);
    
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    Body3dView(const Body3d<c> &sys,
               vector<Body3dState> *xs,
               vector<Vectorcd> *us = 0);
    
    virtual ~Body3dView();
    
  
    virtual void Render(const Body3dState *x,
                        const Vectorcd *u = 0);
    
    void Render(const vector<Body3dState > *xs, 
                const vector<Vectorcd> *us = 0, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Body3d<c> &sys;
    GLUquadricObj *qobj;
    double dirSize;
    
    static void Transform(const Matrix3d &R, const Vector3d &p);
  };
}


using namespace gcop;
using namespace Eigen;

template<int c>
Body3dView<c>::Body3dView(const Body3d<c> &sys) : 
  SystemView<Body3dState, Matrix<double, c, 1> >("Body3d"), sys(sys)
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
Body3dView<c>::Body3dView(const Body3d<c> &sys,
                          vector< Body3dState> *xs,
                          vector<Vectorcd> *us) : 
SystemView<Body3dState, Matrix<double, c, 1> >("Body3d", xs, us), sys(sys)
{
  this->rgba[0] = 0.5;
  this->rgba[1] = 0.5;
  this->rgba[2] = 0.5;
  this->rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
}


template<int c>
Body3dView<c>::~Body3dView()
{
  free(qobj);
}

template<int c>
void Body3dView<c>::Render(const Body3dState *x,
                           const Vectorcd *u)
{
  //   glColor4f(1,0.5,0.5,0.5);
  
  glPushMatrix();
  Transform(x->R, x->p);
  glScaled(sys.ds[0], sys.ds[1], sys.ds[2]); 
  glutSolidCube(1);
  glPopMatrix();
}


template<int c>
void Body3dView<c>::Render(const vector<Body3dState> *xs, 
                           const vector<Vectorcd > *us, 
                           bool rs, 
                           int is, int ie,
                           int dis, int dit,
                           bool dl)
{
  // Viewer::SetColor(this->rgba[0], this->rgba[1], this->rgba[2], this->rgba[3]);
  glColor3f(this->rgba[0], this->rgba[1], this->rgba[2]);

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
    const Body3dState &x = (*xs)[i];
    glVertex3d(x.p[0], x.p[1], x.p[2]);
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
void Body3dView<c>::Transform(const Matrix3d &R, const Vector3d &p)
{
  const SO3 &so3 = SO3::Instance();
  glTranslated(p[0], p[1], p[2]);
  
  Vector3d e;
  so3.log(e, R);
  double n = e.norm();
  if (n > so3.tol) {
    e = e/n;
    glRotated(RAD2DEG(n), e[0], e[1], e[2]);
  }
}


#endif
