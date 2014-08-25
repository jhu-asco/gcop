#include <iomanip>
#include <iostream>
#include "body2dslam.h"
#include "utils.h"
#include "se2.h"
#include "viewer.h"
#include "body2dview.h"
#include "body2dgraphview.h"

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

  Body2d sys(new Body2dForce);
  //  sys.force->D << .01, .01, 1;

  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  Body2dGraph pgt(sys, N, nf);   ///< ground truth
  Body2dGraph pg(sys, N, nf);    ///< noisy pose graph
  pg.odometry = true;

  Body2dGraph::Synthesize3(pgt, pg, tf);

  cout << pg.xs[0].first << endl;
  cout << pg.us[0] << endl;

  Body2dView tview(sys, &pgt.xs);   // true path
  tview.lineWidth = 5;
  tview.rgba[0] = 0;
  tview.rgba[1] = 1;
  tview.rgba[2] = 0;
  Body2dGraphView pgtv(pgt);             // true pose-graph
  pgtv.rgba[0] = 0;   pgtv.rgba[1] = 1;


  Body2dView view(sys, &pg.xs);     // optimized path
  view.lineWidth = 5;
  view.rgba[0] = 1;
  view.rgba[1] = 0;
  view.rgba[2] = 0;

  Body2dGraphView pgv(pg);               // optimized pose-graph
  pgv.rgba[0] = 1; pgv.rgba[0] = 0;

  if (viewer) {
    viewer->Add(tview);
    viewer->Add(pgtv);
    
    viewer->Add(view);
    viewer->Add(pgv);
  }

  struct timeval timer;
  getchar();

  Body2dSlam ba(pg);
  ba.pddp->debug = true; // turn off debug for speed
  ba.pddp->mu = .01;
  ba.pddp->nu = .01;

  for (int i = 0; i < 1000; ++i) {

    cout << "Press Enter to continue" << endl;
    getchar();    
    
    timer_start(timer);
    ba.pddp->Iterate();
    long te = timer_us(timer);
    cout << ba.pddp->dus[0] << endl;
    cout << "Iteration #" << i << " took: " << te << " us." << endl;    
    cout << "p=" << pg.p.transpose() << endl;
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
