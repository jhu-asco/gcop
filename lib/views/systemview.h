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

  using namespace std;

  /**
   * Discrete Dynamical system view.
   */
  template <typename Tx, typename Tu>
    class SystemView : public View {
  public:
    
    /**
     * Create a view for a single state s
     * @param name name
     * @param s rigid body state
     */
    SystemView(const char *name = 0, 
               vector<Tx> *xs = 0,
               vector<Tu> *us = 0);

    virtual ~SystemView();

    virtual void Render(const Tx *x,
                        const Tu *u = 0) = 0;

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
    virtual void Render(const vector<Tx> *xs,
                        const vector<Tu> *us = 0,
                        bool rs = true, 
                        int is = -1, int ie = -1,
                        int dis = 1, int dit = 1,
                        bool dl = true) = 0;

    virtual void Render();
    
    virtual bool RenderFrame(int i);
    
    /**
     * Set the rgb color of the trajectory
     * @param rgb rgb vector
     */
    void SetColor(const double rgba[4]);

    vector<Tx> *xs;  ///< trajectory to render
    vector<Tu> *us;  ///< controls to render
    
    int dis;                   ///< when rendering a trajectory only render a body corresponding to every di-th state (default is 1)
    
    int dit;                   ///< when rendering a trajectory only put a point at every di-th state (default is 1)
    
    double rgba[4];            ///< color

    double lineWidth;         ///< trajectory line width (default is 1)

    bool renderSystem;        ///< whether to render the physical system (e.g. a body)
    bool renderForces;       ///< whether to render forces
  };

  
  template <typename Tx, typename Tu>
    SystemView<Tx,Tu>::SystemView(const char *name, vector<Tx> *xs, vector<Tu> *us) : View(name),
    xs(xs),
    us(us),
    dis(1),
    dit(1),
    lineWidth(1),
    renderSystem(true),
    renderForces(false) {
    rgba[0] = 1;
    rgba[1] = 1;
    rgba[2] = 1;
    rgba[3] = 0;
  }
  
  template <typename Tx, typename Tu>
    SystemView<Tx, Tu>::~SystemView() {
  }
  
  template <typename Tx, typename Tu>
    void SystemView<Tx,Tu>::SetColor(const double rgba[4]) {
    memcpy(this->rgba, rgba, 4*sizeof(double));
  }
  
  template <typename Tx, typename Tu>
    void SystemView<Tx, Tu>::Render() {
    glColor4dv(rgba);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //  glEnable(GL_BLEND);
    
    //  glColor4f(1,0,0,0);
    
    if (xs)
      Render(xs, us, renderSystem, 0, xs->size()-1, dis, dit); 
    
    //  glDisable(GL_BLEND);
  }
  
  
  template <typename Tx, typename Tu>
    bool SystemView<Tx, Tu>::RenderFrame(int i) {
    
    // glColor4dv(rgba);
    
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //  glEnable(GL_BLEND);
    
    if (xs) {
      if (i < 0 || i >= xs->size())
        return false;
      
      // draw path without states (only last if applicable)
      Render(xs, us, false, 0, i, dis, dit, false);
      
      // draw current state
      if (us && i < us->size())
        Render(&(*xs)[i], &(*us)[i]);
      else
        Render(&(*xs)[i]);
        
      //    glDisable(GL_BLEND);
      // more states left
      return true;
    }
    
    //  glDisable(GL_BLEND);
    
    return false;
  }
}

#endif
