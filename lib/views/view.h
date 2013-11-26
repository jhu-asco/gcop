#ifndef GCOP_VIEW_H
#define GCOP_VIEW_H

#include <pthread.h>

namespace gcop {

class Viewer;

/**
 * Basic view class that all views should subclass. Render() should provide the rendering. 
 * Lock() and Unlock() should be used when any data used in Render() is being updated.
 *
 * Author: Marin Kobilarov -- Copyright (C) 2005
 */
class View {
 public:
  /**
   * Create a view. All subclasses should call the constructor
   * @param name optional name
   */
  View(const char *name = 0);

  /**
   * Destroy view
   */
  virtual ~View();

  /**
   * Render a view: all subclasses should provide this method
   */
  virtual void Render();

  /**
   * Render a view: all subclasses should provide this method
   * @param i frame #
   * @return return true if more frames should be rendered
   */
  virtual bool RenderFrame(int i = 0);  

  /**
   * Lock a view. This is done when updating data inside this object. 
   * While a view is locked Render() is not called
   */
  void Lock();

  /**
   * Unlock the view and continue display
   */
  void Unlock();

  /**
   * The view can optionally specify 
   */
  //  virtual bool GetCamera(double x[3], double rpy[3]) const;
  
 protected:
  
  char name[256];         ///< name
  int id;                 ///< id
  pthread_mutex_t mut;    ///< synchronization

  friend class Viewer;
  
}; 

}

#endif
