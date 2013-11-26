#include "ce.h"
#include <iostream>
#include "utils.h"
#include <assert.h>

using namespace Eigen;
using namespace gcop;
using namespace std;

// produce a discrete trajectory xs from parameter z, starting at initial state x0
void Traj(vector<Vector2d> &xs, const VectorXd &z, const Vector2d &x0)
{
  int K = z.size()/2; // trajectory segments
  assert(xs.size() == K+1);

  xs[0] = x0;
  for (int k = 0; k < K; ++k) {
    const Vector2d &v = z.segment(2*k,2);  // velocity
    xs[k+1] = xs[k] + v;                   // current state
  }
}


// compute the cost of a parameter z, with initial state x0 and goal xf
double Cost(const VectorXd &z, const Vector2d &x0, const Vector2d &xf)
{
  int K = z.size()/2;   // trajectory segments
  
  Vector2d Qf(1, 1);     // LQR-type final cost matrix (diagonal)
  Vector2d R(.01, .01);  // LQR-type control cost matrix (diagonal)
  
  double c = 0;
  Vector2d x = x0;

  for (int k = 0; k < K; ++k) {
    const Vector2d &v = z.segment(2*k,2);  // velocity
    x += v;                                // current state
    c += v.dot(R.cwiseProduct(v))/2;       // total distance
  }
  Vector2d dx = x - xf;                    // vector to goal
  c += dx.dot(Qf.cwiseProduct(dx))/2;      // final cost
  
  return c;
}



int main(int argc, char** argv)
{
  int K = 5;
  
  MatrixXd S = VectorXd::Constant(2*K, .01).asDiagonal();

  Ce ce(2*K, 1, &S);

  // initialize using unit variance centered at zero
  ce.gmm.ns[0].mu.setZero();
  ce.gmm.ns[0].P = VectorXd::Constant(2*K, 1).asDiagonal();
  ce.gmm.Update();

  // initial and final positions
  Vector2d x0(0,0);
  Vector2d xf(5,5);

  struct timeval timer;

  // number of samples
  int N = 100;
  
  VectorXd z(2*K);            // parameter vector
  vector<Vector2d> xs(K+1);   // corresponding trajectory

  while(1) {
    timer_start(timer);
    ce.Reset();
    
    // add samples
    for (int j = 0; j < N; ++j) {
      ce.Sample(z);
      ce.AddSample(z, Cost(z, x0, xf));
    }

    // estimate distribution
    ce.Select();    
    if (!ce.Fit()) {
      cout << "[W] TrajectoryPrmSample::Sample: ce.Fit failed!" << endl;
      break;
    }
   
    long ct = timer_us(timer);

    // construct trajectory using the first sample (this is the one with lowest cost)
    Traj(xs, ce.zps[0].first, x0);
    cout << "Solution: xs[K]=" << xs[K].transpose() << " c=" << ce.cs[0] << endl;

    cout << "Took: " << ct << " us. Press Enter to continue..." << endl;
    getchar();

  }
  return 0;
}  
