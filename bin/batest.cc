#include <iomanip>
#include <iostream>
#include "ba.h"
#include "viewer.h"
#include "kinbody2dview.h"
#include "utils.h"
#include "se2.h"
#include "posegraph2dview.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

void Run(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(18.875, 1.625, -4, -1.75, -5.45);

  int N = 50;
  double tf = 25;
  double h = tf/N;
  int nf = 5*(N+1);

  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  Posegraph2d pgt(N, nf);   ///< ground truth
  Posegraph2d pg(N, nf);    ///< noisy pose graph
  Posegraph2d::Synthesize2(pgt, pg, tf);

  Ba ba(pg);

  Kinbody2dView tview(ba.sys, &pgt.gs);   // true path
  tview.lineWidth = 5;
  tview.rgba[0] = 0;
  tview.rgba[1] = 1;
  tview.rgba[2] = 0;
  Posegraph2dView pgtv(pgt);             // true pose-graph
  pgtv.rgba[0] = 0;   pgtv.rgba[1] = 1;

  Kinbody2dView view(ba.sys, &pg.gs);     // optimized path
  view.lineWidth = 5;
  view.rgba[0] = 1;
  view.rgba[1] = 0;
  view.rgba[2] = 0;

  Posegraph2dView pgv(pg);               // optimized pose-graph
  pgv.rgba[0] = 1; pgv.rgba[0] = 0;


  if (viewer) {
    viewer->Add(tview);
    viewer->Add(pgtv);
    
    viewer->Add(view);
    viewer->Add(pgv);
  }

  struct timeval timer;

  //  ba.pddp.debug = false; // turn off debug for speed
  ba.pddp.nu = .1;

  for (int i = 0; i < 1000; ++i) {

    cout << "Press Enter to continue" << endl;
    getchar();    
    
    timer_start(timer);
    ba.pddp.Iterate();
    long te = timer_us(timer);
    
    cout << "Iteration #" << i << " took: " << te << " us." << endl;    
  }
  
  cout << "done!" << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{

#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "videos/sys";

  pthread_t dummy;
  pthread_create( &dummy, NULL, (void *(*) (void *)) Run, viewer);

#else
  Run(0);
#endif


#ifdef DISP
  viewer->Start();
#endif


  return 0;
}
