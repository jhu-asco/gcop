#include "imuukf.h"
#include <iostream>
#include "quat.h"
#include <limits>
#include <sys/time.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace gcop;
using namespace Eigen;

/**
 * 0-mean unit-stdev normal distribution 1-D sample
 * @return sample
 */
double randn()
{
  static const double r_max=RAND_MAX;
  double U1,U2,V1,V2;
  double S = 2;
  while(S>=1) {
    U1 = rand()/r_max;
    U2 = rand()/r_max;
    V1 = 2*U1-1;
    V2 = 2*U2-1;
    S = V1*V1+V2*V2;
  }
  return V1*sqrt(-2*log(S)/S);
}


#define RND (rand()/(double)RAND_MAX - .5)

/**
 * Start a timer
 * @param timer timer
 */
inline void timer_start(struct timeval* timer)
{
  gettimeofday(timer, 0);
}

/**
 * Get elapsed time in microseconds
 * Timer should be started with timer_start(timer)
 * @param timer timer
 * @return elapsed time
 */
inline long timer_us(struct timeval* timer)
{
  struct timeval now;
  gettimeofday(&now, 0);
  return (now.tv_sec - timer->tv_sec)*1000000 + now.tv_usec - timer->tv_usec;
}


int main(int argc, char** argv)
{
  ImuModel im;
  ImuUKF ukf(im);
  int N = 50;
  
  Vector3d wt;
  wt << .2, .4, .1;

  VectorXd xt(10);
  xt.setZero();
  xt[0] = 1;
  VectorXd ut(6);           
  
  double dt = .1;
  im.dt = dt;

  struct timeval timer;

  for (int i = 0; i < N; ++i) {
    double t = i*dt;
    Quat qi(xt.data());
    qi.Invert();

    // observed
    Vector3d w(3);
    Vector3d a(3);
    Vector3d m(3);
    
    // true 
    Vector3d at(3);
    Vector3d mt(3);
    qi.Rotate2(at.data(), im.g.data());
    qi.Rotate2(mt.data(), im.m0.data());

    // generate noisy
    for (int j = 0; j < 3; ++j) {
      w[j] = wt[j] + im.wn[j]*randn();
      a[j] = at[j] + im.vn[j]*randn();
      m[j] = mt[j] + im.mon[j]*randn();
    }
    timer_start(&timer);

    cout << w << " " << a << " " << m << endl;
    ukf.Process(t, w, a, m);

    long rt = timer_us(&timer);
    cout << "time=" << rt << " us." << endl;
    

    if (i > 0) {
      memcpy(ut.data(), wt.data(), 3*sizeof(double));
      memcpy(ut.data() + 3, at.data(), 3*sizeof(double));
      
      VectorXd xn(xt.size());
      im.f(xn, xt, ut);
      xt = xn;
      
      cout << "True: " << xt << endl;
      cout << "Estm: " << ukf.x << endl;
      
      cout << "L_2 norm (state) =" << (xt - ukf.x).norm() << endl;
      cout << "L_2 norm (quat) =" << (xt.head(4) - ukf.x.head(4)).norm() << endl;
    }
  }
}
