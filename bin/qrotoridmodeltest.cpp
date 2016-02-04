#include <iomanip>
#include <iostream>
#include "utils.h"
#include "qrotoridmodelcost.h"
#include "params.h"


//#define USE_SDDP

#ifdef USE_SDDP
#include "sddp.h"
#else
#include "ddp.h"
#endif

using namespace std;
using namespace Eigen;
using namespace gcop;


#ifdef USE_SDDP
typedef SDdp<QRotorIDState, 15, 4, 10> QrotorDdp;
#else
typedef Ddp<QRotorIDState, 15, 4, 10> QrotorDdp;
#endif

Params params;

void solver_process()
{

  int N = 32;      // discrete trajectory segments
  double tf = 2;   // time-horizon

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);  

  double h = tf/N; // time-step

  // system
  QRotorIDModel sys;


  QRotorIDState x0;
  VectorXd qv0(15);
  params.GetVectorXd("x0", qv0);  
  SO3 &so3 = SO3::Instance();
  so3.q2g(x0.R, qv0.head(3));
  x0.p = qv0.segment<3>(3); x0.w = qv0.segment<3>(6); x0.v = qv0.segment<3>(9); x0.u =qv0.tail<3>();

  QRotorIDState xf;
  VectorXd qvf(12);
  params.GetVectorXd("xf", qvf);  
  so3.q2g(xf.R, qvf.head(3));
  xf.p = qvf.segment<3>(3); xf.w = qvf.segment<3>(6); xf.v = qvf.tail<3>(); xf.u.setZero();

  QRotorIDModelCost cost(sys,tf,xf);

  VectorXd Q(15);
  VectorXd R(4);  
  VectorXd Qf(15);

  params.GetVectorXd("Q", Q);  
  params.GetVectorXd("R", R);
  params.GetVectorXd("Qf", Qf);

  cost.Q = Q.asDiagonal();
  cost.R = R.asDiagonal();
  cost.Qf = Qf.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // states
  vector<QRotorIDState> xs(N+1, x0);

  // initial controls (e.g. hover at one place)
  vector<Vector4d> us(N);
  for (int i = 0; i < N; ++i) {
    us[i].tail(3).setZero();
    us[i][0] = (9.81/sys.kt);
    us[i][3] =  0.05*cos((double(i)/N)*2*M_PI);
  }
    
  QrotorDdp ddp(sys, cost, ts, xs, us);
  ddp.mu = .01;

  /*struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed
  for (int i = 0; i < 100; ++i) {
    timer_start(timer);
    ddp.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << ": took " << te << " us." << endl;        
    //getchar();
  }
  */

  //Display Current states:
  cout<<"Index Pos Vel rpy Omega commandedrpy"<<endl;
  for(int i = 0; i < N; i++)
  {
    //  sys.Step(xs[i+1],i*h,xs[i],us[i],h,0,0,0,0);
      Vector3d rpy;
      so3.g2q(rpy, xs[i].R);
      cout<<i<<" "<<xs[i].p.transpose()<<" "<<xs[i].v.transpose()<<" "<<rpy.transpose()<<" "<<xs[i].w.transpose()<<" "<<xs[i].u.transpose()<<" "<<us[i].transpose()<<endl;
  }

  cout << "done!" << endl;
  cout << "Click on window and press 'a' to animate." << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{
  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("qrotormodelid.cfg");


  solver_process();


  return 0;
}
