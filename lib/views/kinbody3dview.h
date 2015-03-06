#ifndef GCOP_KINBODY3DVIEW_H
#define GCOP_KINBODY3DVIEW_H

#include "kinbody3d.h"
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
  
    template <int _nu = 6>
    class Kinbody3dView : public SystemView<Matrix4d, Matrix<double, _nu, 1> > {

    typedef Matrix<double, 6, 1> Vector6d;
    typedef Matrix<double, _nu, 1> Vectorud;

  public:
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     */
    
    Kinbody3dView(const Kinbody3d<_nu> &sys);
    
    
    /**
     *  Create a particle view of trajectory traj
     * @param sys particle
     * @param xs trajectory
     */
    Kinbody3dView(const Kinbody3d<_nu> &sys,
               vector<Matrix4d > *xs,
               vector<Vectorud> *us = 0);
    
    virtual ~Kinbody3dView();
    
  
    virtual void Render(const Matrix4d *x,
                        const Vectorud *u = 0);
    
    void Render(const vector<Matrix4d > *xs, 
                const vector<Vectorud> *us = 0, 
                bool rs = true,
                int is = -1, int ie = -1,
                int dis = 1, int dit = 1,
                bool dl = false);
    
    const Kinbody3d<_nu> &sys;
    GLUquadricObj *qobj;
    double dirSize;
    
    static void Transform(const Matrix3d &R, const Vector3d &p);
  };
}


using namespace gcop;
using namespace Eigen;

template <int _nu>
Kinbody3dView<_nu>::Kinbody3dView(const Kinbody3d<_nu> &sys) : 
  SystemView<Matrix4d, Matrix<double, _nu, 1> >("Kinbody3d"), sys(sys)
{
  this->rgba[0] = .5;
  this->rgba[1] = .5;
  this->rgba[2] = .5;
  this->rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
  dirSize = -1;
}

template <int _nu>
Kinbody3dView<_nu>::Kinbody3dView(const Kinbody3d<_nu> &sys,
                          vector< Matrix4d > *xs,
                          vector<Vectorud> *us) : 
SystemView<Matrix4d, Matrix<double, _nu, 1> >("Kinbody3d", xs, us), sys(sys)
{
  this->rgba[0] = 0.5;
  this->rgba[1] = 0.5;
  this->rgba[2] = 0.5;
  this->rgba[3] = 0;
  this->lineWidth = 2;
  qobj = gluNewQuadric();
}


template <int _nu>
Kinbody3dView<_nu>::~Kinbody3dView()
{
  free(qobj);
}

template <int _nu>
void Kinbody3dView<_nu>::Render(const Matrix4d *x,
                           const Vectorud *u)
{
  //   glColor4f(1,0.5,0.5,0.5);
  
  glPushMatrix();
  Transform(x->block<3,3>(0,0), x->block<3,1>(0,3));
  glScaled(sys.d[0], sys.d[1], sys.d[2]); 
  glutSolidCube(1);
  glPopMatrix();
}


template <int _nu>
void Kinbody3dView<_nu>::Render(const vector<Matrix4d > *xs, 
                           const vector<Vectorud > *us, 
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
    const Matrix4d &x = (*xs)[i];
    glVertex3d(x(0,3), x(1,3), x(2,3));
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

template <int _nu>
void Kinbody3dView<_nu>::Transform(const Matrix3d &R, const Vector3d &p)
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
