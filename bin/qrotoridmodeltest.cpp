#include <iomanip>
#include <iostream>
#include "utils.h"
#include "qrotoridgndocp.h"
#include "qrotoridmodelcost.h"
#include "params.h"
#include "unistd.h"


#define USE_SPLINEPARAM

#ifdef USE_SPLINEPARAM
#include "uniformsplinetparam.h"
#else
#include "controltparam.h"
#endif


using namespace std;
using namespace Eigen;
using namespace gcop;



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
  VectorXd x0_std(15);
  VectorXd p_std(13);
  VectorXd obs(7);
  Matrix<double,13,1> p_mean;

  params.GetVectorXd("Q", Q);  
  params.GetVectorXd("R", R);
  params.GetVectorXd("Qf", Qf);
  params.GetVectorXd("obs",obs);

  cost.Q = Q.asDiagonal();
  cost.R = R.asDiagonal();
  cost.Qf = Qf.asDiagonal();

  vector<Obstacle> obstacles(1);
  obstacles[0].set(obs);

  cost.UpdateGains();

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
    us[i][0] = (9.81/sys.kt)+0.1;
  }

  int Nk = N/5;//Number of segments
  VectorXd tks(Nk+1);
  //vector<double> tks(Nk+1);
  for (int k = 0; k <=Nk; ++k)
    tks[k] = k*(tf/Nk);

#ifdef USE_SPLINEPARAM
  UniformSplineTparam<QRotorIDState, 15, 4,13 > ctp(sys, tks,2);//Second order spline
#else
  ControlTparam<QRotorIDState, 15, 4, 13> ctp(sys, tks);
#endif

  p_mean<<sys.kt, sys.kp, sys.kd, sys.a0, sys.tau0;//Straight params from system

  QRotorIdGnDocp gn(sys, cost, ctp, ts, xs, us, &p_mean);

  /*for(int i = 0; i < (N+1); i++)
  {
    cout<<ts[i]<<"\t"<<xs[i].p.transpose()<<endl;
  }
  return;
  */


  //Generate Stdev of Trajectories using current control:
  params.GetVectorXd("x0_std", x0_std);
  params.GetVectorXd("p_std", p_std);
  gn.stdev_initial_state = x0_std.asDiagonal();
  gn.stdev_params = p_std.asDiagonal();

  //Add one obstacle:
  gn.AddObstacles(obstacles);

  gn.GenerateStdev();

  cout<<"Generated Stdev: "<<endl;
  /*for(int i = 0; i <=N; i++)
  {
     cout<<ts[i]<<"\t"<<gn.sample_mean_params.col(i).transpose()<<"\t"<<gn.xs_std.col(i).transpose()<<endl;
  }

  return;////DEBUG
  */


  struct timeval timer;
  //  ddp.debug = false; // turn off debug for speed
  for (int i = 0; i < 10; ++i) {
    timer_start(timer);
    gn.Iterate();
    gn.ko = 5*gn.ko;
    gn.Reset();
    long te = timer_us(timer);
    cout<<gn.J<<endl;
//    cout << "Iteration #" << i << ": took " << te << " us." << endl;        
    //getchar();
  }

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
