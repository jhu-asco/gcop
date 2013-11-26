#include <iostream>
#include "viewer.h"
#include "utils.h"
#include <pthread.h>
#include <assert.h>
#include "gmmview.h"

using namespace dgc;
using namespace std;
using namespace itpp;


void solver_process(Viewer* viewer)
{
  Gmm gmm(2,2);
  gmm.Init(vec("0.0, 0.0"), vec("6.0 ,6.0"));

  Normal n1(vec("1.0,1"), diag(vec("2.0,1")));
  Normal n2(vec("5.0,5"), diag(vec("1.0,2")));

  int N = 3;
  vec* xs[N];
  for (int j = 0; j < N; ++j) {
    xs[j] = new vec(2);
    if (RND < .5)
      n1.Sample(*xs[j]);
    else
      n2.Sample(*xs[j]);
  }

  *xs[0] = vec("1 1");
  //  *xs[1] = vec("1.3 1.2");
  *xs[2] = vec("5 5.1");
  *xs[1] = vec("4.5 5");

  GmmView gv(gmm, 2);
  viewer->AddView(gv);
  
  mat S = diag(vec("0.05 0.05"));

  getchar();  

  
  struct timeval timer;
  for (int i = 0; i >= 0; ++i) {
    gv.Lock();        
    timer_start(&timer);
    gmm.Fit((const vec**)xs, N, 0, 1, &S);
    cout << i << ": took " << timer_us(&timer) << " us" << endl;
    gv.Unlock();
    
    getchar();
  }
}


#define DISP


int main(int argc, char** argv)
{

#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "pnt/frames/pnt";
  viewer->displayName = "pnt/display/pnt";

  pthread_t dummy;
  pthread_create( &dummy, NULL, (void *(*) (void *)) solver_process, viewer);

#else
  solver_process(0);
#endif


#ifdef DISP
  viewer->Start();
#endif


  return 0;
}  
