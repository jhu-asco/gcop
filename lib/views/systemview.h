#ifndef GCOP_SYSTEMVIEW_H
#define GCOP_SYSTEMVIEW_H

#include "GL/glu.h"
#include "GL/glut.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "viewer.h"
#include "utils.h"
#include "system.h"
#include "view.h"


namespace gcop {

  /**
   * Discrete Dynamical system view.
   */
  template <typename T>
    class SystemView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    SystemView(const char *name = 0, vector<T> *xs = 0);

    virtual ~SystemView();

    virtual void Render(const T &x) = 0;

    /**
     * Render a trajectory
     * @param traj trajectory
     * @param rs render state
     * @param is start index
     * @param ie end index
     * @param dis index step for rendering bodies along trajectory
     * @param dit index step for rendering states(points) along trajectory
     * @param dl always draw last state
     */
    virtual void Render(const vector<T>& xs,
                        bool rs = true, 
                        int is = -1, int ie = -1,
                        int dis = 1, int dit = 1,
                        bool dl = true) = 0;

    /**
     * Render only states along the trajectory at time-steps h
     * @param tarj trajectory
     * @param h time-step
     */
    void RenderStates(const vector<T>& xs, double h);    

    virtual void Render();
    
    virtual bool RenderFrame(int i);
    
    /**
     * Set the rgb color of the trajectory
     * @param rgb rgb vector
     */
    void SetColor(const double rgba[4]);

    vector<T> *xs;  ///< trajectory to render
    
    int dis;                   ///< when rendering a trajectory only render a body corresponding to every di-th state (default is 1)
    
    int dit;                   ///< when rendering a trajectory only put a point at every di-th state (default is 1)
    
    double rgba[4];            ///< color

    double lineWidth;         ///< trajectory line width (default is 1)
  };

  
  template <typename T>
    SystemView<T>::SystemView(const char *name, vector<T> *xs) : View(name),
    xs(xs),
    dis(1),
    dit(1),
    lineWidth(1) {
    rgba[0] = 1;
    rgba[1] = 1;
    rgba[2] = 1;
    rgba[3] = 0;
  }
  
  template <typename T>
    SystemView<T>::~SystemView() {
  }
  
  template <typename T>
    void SystemView<T>::SetColor(const double rgba[4]) {
    memcpy(this->rgba, rgba, 4*sizeof(double));
  }
  
  

  
  template <typename T>
    void SystemView<T>::RenderStates(const vector<T> &xs, double h) {
    /*
      assert(h>0);
      glColor4dv(rgba);
      State s(traj.sys);
      for (s.t =0 ; s.t <= traj.states[traj.sn]->t; s.t+=h) {
      traj.Get(s);
      Render(s);  
      }
    */
  }
  
  
  template <typename T>
    void SystemView<T>::Render() {
    glColor4dv(rgba);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //  glEnable(GL_BLEND);
    
    //  glColor4f(1,0,0,0);
    
    if (xs)
      Render(*xs, true, 0, xs->size()-1, dis, dit); 
    
    //  glDisable(GL_BLEND);
  }
  
  
  template <typename T>
    bool SystemView<T>::RenderFrame(int i) {
    
    // glColor4dv(rgba);
    
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //  glEnable(GL_BLEND);
    
    if (xs) {
      if (i < 0 || i >= xs->size())
        return false;
      
      // draw path without states (only last if applicable)
      Render(*xs, false, 0, i, dis, dit, false);
      
      // draw current state
      Render((*xs)[i]);
      
      //    glDisable(GL_BLEND);
      // more states left
      return true;
    }
    
    //  glDisable(GL_BLEND);
    
    return false;
  }
}

#endif
