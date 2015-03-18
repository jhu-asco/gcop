#ifndef GCOP_VIEWER_H
#define GCOP_VIEWER_H

#include <stdlib.h>
#include <GL/glut.h>
#include <map>
#include "view.h"

namespace gcop {


/**
 *  Light-weight OpenGL viewer
 *  The viewer keeps a list of "views" of different
 *  objects and displays them.
 *  A view is simply defined as an object that "knows"
 *  how to be rendered: it provides its own Render() function
 *
 *  Author: Marin Kobilarov (mkobilar(at)robotics.usc.edu)
 */
class Viewer {
 public:
  Viewer();
  virtual ~Viewer();

  /**
   *  Initialize the viewer with the command line args
   */
  void Init(int* argc, char** argv);

  /**
   *  Start the viewer and loop until the program is killed
   *  NOTE: This should be done in the main process thread and not in a seperate thread!
   */
  void Start();

  /**
   *  Add a view
   */
  void Add(View& view);

  /**
   *  Remove a view
   */
  void Remove(View& view);


  /**
   *  Render all views. This is automatically called by graphics engine.
   *  You can call this method only if you want to render the views immediately
   */
  void RenderViews();

  /**
   *  Draw origin. 
   */
  void DrawOrigin(GLUquadricObj *qobj = 0);

  /**
   * Draw frame at the origin with unit length axes
   */
  static void DrawFrame(GLUquadricObj *qobj = 0);

  /** 
   * Draws an arrow with length l
   * @param l length of arrow
   * @param qobj quadratic object used for drawing (leave 0 to create a new one)
   */
  static void DrawArrow(double l = 1.0, GLUquadricObj *qobj = 0);

  /** 
   * Draws an arrow aligned with vector v
   * @param v 3x1 vector
   * @param qobj quadratic object used for drawing (leave 0 to create a new one)
   */
  static void DrawArrow(const double v[3], GLUquadricObj *qobj = 0);


  /**
   *  Draws text
   *  @param str string
   *  @param size font size
   */
  static void DrawText(const char *str, float size = .2);

  /**
   *  Draw all views and origin  
   */
  void DrawAll();

  /**
   * Set the default background (clear) color
   * @param rgb rgb color
   */
  void SetColor(const float rgb[3]);

  float sphi;
  float stheta;
  float sx, sy, sz;
  float zNear;
  float zFar;
  float aspect;
  float xcam;
  float ycam;
  int width;
  int height;

  float os; ///< origin scale (how long are the origin axes, default is 1)

  float rgb[3]; ///< default background (clear) color

  void Reshape(int width, int height); 
  void Render();
  void Keyboard(unsigned char ch, int x, int y);
  void Mouse(int button, int state, int x, int y);
  void Motion(int x, int y);

  /**
   * Save an image of the current display
   * @param name file name (if omited then the name will be
   *         "%s_%04d.ppm" where %s is Viewer::frameName and i
   *         is a global index
   */
  void SaveImage(const char* name = 0);

  static GLuint LoadTexture(const char * fname, int wrap = true);

  /**
   * Draw a circle
   * @param x x-coord
   * @param y y-coord
   * @param r radius
   * @param slices resolution in slices
   */
  static void DrawCircle(double x = 0, double y = 0, double r = 1, int slices = 20);
  
  int window;
 protected:

  int downX;
  int downY;
  int leftButton;
  int middleButton;
  int rightButton;

  pthread_mutex_t mut;
  
  std::map<int, View*> views;

  GLUquadricObj *aqobj;

 public:
  bool animate;                 ///< whether to animate all trajectory frames
  double animDelay;             ///< animation delay (def. is 30 ms)

  bool saveFrames;              ///< whether to save frame images
  const char* frameName;        ///< frame image filename

  bool saveDisplay;             ///< whether to continuously save display images at every refresh
  const char* displayName;      ///< display image filename
  int displayNo;                ///< counter for display image filenames

  bool saveSnapshot;            ///< save a single snapshot

  bool altLight;

  /**
   * Draw a cylinder between two points
   * @param xa starting point
   * @param xb end point
   * @param r radius
   * @param sd unused
   */
  static void Cylinder(const double xa[3], const double xb[3], 
                       double r, int sd);
    
  void SetCamera(float sphi, float stheta,
                 float sx, float sy, float sz);

  void SetCamera(const float params[5]);


  static void SetColor (float r, float g, float b, float alpha = 0);  

  static void SetMaterial( GLfloat ambientR, GLfloat ambientG, GLfloat ambientB, 
                           GLfloat diffuseR, GLfloat diffuseG, GLfloat diffuseB, 
                           GLfloat specularR, GLfloat specularG, GLfloat specularB,
                           GLfloat shininess );
    
  static Viewer* instance;

 protected:
  unsigned char *rgbimg;    ///< rgb image used internally for saving
  unsigned char *ppmimg;    ///< ppm image used internally for saving  
};

}
#endif


/*! \mainpage OpenGL Viewer
 *  
 * \section Intro 
 *
 *  OpenGL viewer
 *  The viewer keeps a list of "views" of different
 *  objects and displays them.
 *  A view is simply defined as an object that "knows"
 *  how to be rendered: it provides its own Render() function
 *
 * \section Installation
 * \subsection Requirements
 *  Linux/Unix, g++ compiler, opengl
 *
 * \subsection Download
 *  http://www.cds.caltech.edu/~marin/projects/viewer/viewer-0.1.tar.gz
 *
 * \subsection Instructions
 * - download: viewer-0.1.tar.gz
 * - unpack: tar xfz viewer-0.1.tar.gz
 * - type: cd viewer-0.1, make install
 * - to test: cd test, ./viewertest
 *
 * \section Usage
 * see test executable
 *
 * \section Author
 *  Marin Kobilarov -- Copyright (C) 2005
 * \section Keywords
 * opengl, viewer
 *
 */


