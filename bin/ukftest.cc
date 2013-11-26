#include "ukf.h"
#include <iostream>

using namespace std;
using namespace gcop;
using namespace Eigen;

class ParticleModel : public Model {
public:
  ParticleModel() :
    Model(2, 1), dt(0.1) {
    double s = .5;
    Vector2d w(2);
    w(0) = dt*dt/2*s;
    w(1) = dt*s;
    
    Q.diagonal() << w.cwiseProduct(w);
    
    R(0,0) = 1;
  }
    
  bool f (VectorXd &xn,
          const VectorXd &x)
  {
    xn[0] = x[0] + dt*x[1];
    xn[1] = x[1];
    return true;
  }
  
  bool h(VectorXd &z, 
         const VectorXd &x)
  {
    z[0] = x[0];
    return true;
  }  
  
  double dt; /// < time-step
};


int main(int argc, char** argv)
{
  ParticleModel pm;
  UKF ukf(pm);

  int N = 100;
  
  VectorXd xn(2);
  VectorXd xt(2);
  xt << 0, 1;

  ukf.x = xt;
  ukf.P = pm.Q;

  cout << pm.Q << endl;
  cout << ukf.P << endl;
  cout << ukf.x << endl;

  VectorXd z(1);

  for (int i = 0; i < N; ++i) {
    ukf.Predict();
    
    pm.f(xn, xt);
    xt = xn;

    z[0] = xt[0] + sqrt(pm.R(0,0))*(rand()/(double)RAND_MAX - .5);
    
    ukf.Update(z);
    
    cout << "True: " << xt[0] << "," << xt[1] << endl;
    cout << "Estm: " << ukf.x[0] << "," << ukf.x[1] << endl;
  }
}
