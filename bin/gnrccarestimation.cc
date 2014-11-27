#include <iomanip>
#include <iostream>
//#include "gndocp.h"
#include "viewer.h"
#include "rccarview.h"
#include "utils.h"
#include "lqsensorcost.h"
#include "params.h"
#include "point3dgps.h"
#include "gndoep.h"

//#include "controltparam.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

//typedef  RccarGnest;
//typedef Matrix<double, 1, 1> Vector1d;
//typedef Matrix<double, 9, 1> Vector9d;

Params params;

void projectmanifold(const Vector4d &rccarstate, Point3dState &pstate)
{
  pstate.q.head<2>() = rccarstate.head<2>();//First 2 elements are x and y
  pstate.q[2] = 0;//Set z to 0. This is enough for GPS
}
/*
inline void projectmanifold(const Point3dState &pstate, Vector4d &rccarstate)
{
  rccarstate.head<2>() = pstate.q.head<2>(); //Does not correct the angle of the car or steering angle
}
*/

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(-5, 51, -0.2, -0.15, -2.3);

  int N = 3;        // number of segments

  int iters = 3;

  params.GetInt("N", N);  

  params.GetInt("iters", iters);
  




  //  sys.U.lb[1] = tan(-M_PI/5);
  //  sys.U.ub[1] = tan(M_PI/5);

  // initial state
  Vector4d x0(1,1,0,0);
  params.GetVector4d("x0", x0);

  // final state
  //Vector4d xf(0,0,0,0);
  //params.GetVector4d("xf", xf);  

  Rccar sys;
  //Add sensor
  Point3dGps<2> gps;//Gps sensor with controls and nofparams
  // cost
  LqSensorCost<Vector4d, 4, 2, Dynamic, 9, Vector3d, 3> cost(sys, gps.Z);

  VectorXd R(3);
  if (params.GetVectorXd("R", R))
    cost.R = R.asDiagonal();

  VectorXd S(4);
  if (params.GetVectorXd("S", S))
    cost.S = S.asDiagonal();
  
  VectorXd P(2);
  if (params.GetVectorXd("P", P)) 
    cost.P = P.asDiagonal();

  cost.UpdateGains();

  // times
  vector<double> ts(N+1);
  // states
  vector<Vector4d> xs(N+1);
  vector<Vector3d> zs(N+1);
  vector<Vector2d> us(N);

  VectorXd mup(2);//Initial Prior
  VectorXd p0(2);//Initial Guess

  //mup[0] = 0.4;//True value is 0.3 length of car testing

  //p0[0] = 1.6;//Initial guess for parameter
  if (!params.GetVectorXd("p0", p0)) 
  {
    cout<<"Cannot find p0 initial guess for parameters"<<endl;
    return;
  }
  mup = p0;//Copy the initial guess to be the same as prior for the parameters
 
  //Check if zs, us, ts are provided
  VectorXd ts1(N+1);
  if(params.GetVectorXd("ts",ts1))
  {
    VectorXd zs1(2*(N+1));
    //zs1<<0,0.0059947,0.023879,0.053635,0.095283,0.14884,0.21425,0.29147,0.38039,0.48086,0.59293,0.71653,0.8509,0.99614,1.1513,1.3165,1.4904,1.673,1.8632,2.0603,2.2639,2.473,2.6871,2.9055,3.1275,3.3526,3.5812,3.815,4.056,4.3054,4.5648,4.8355,5.1184,5.4146,5.7245,6.0484,6.3861,6.7367,7.0986,7.4694,7.8462,0,6.2752e-07,1.5211e-06,1.706e-05,9.0831e-05,0.00032536,0.0009157,0.0021838,0.0046084,0.0088099,0.015545,0.025724,0.04043,0.06081,0.087817,0.12265,0.16556,0.21787,0.28003,0.35252,0.43639,0.5325,0.63997,0.75866,0.88852,1.0291,1.1787,1.3338,1.4912,1.6485,1.803,1.9521,2.0934,2.2238,2.3402,2.4394,2.5175,2.5707,2.595,2.5864,2.5408;
    //[Guess 0.62, 1.602]


    VectorXd us1(2*N);
    //us1<<0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0,0.0012039,0.0048331,0.0099351,0.015952,0.023033,0.030793,0.038543,0.045336,0.051973,0.057203,0.061823,0.062149,0.063442,0.060825,0.064143,0.062277,0.067747,0.071337,0.075629,0.076314,0.07029,0.056808,0.038404,0.013267,-0.01683,-0.042539,-0.058752,-0.07284,-0.084109,-0.093927,-0.10311,-0.11189,-0.12043,-0.12874,-0.13685,-0.14475,-0.15251,-0.16001,-0.16716;

    ts[0] = ts1[0];
    for(int i = 0;i < N;i++)
    {
      zs[i]<<zs1[i], zs1[i+N+1], 0;
      us[i]<<us1[i], tan(us1[i+N]);
      ts[i+1] = ts1[i+1];
      cout<<"zs["<<i<<"]: "<<zs[i].transpose()<<endl;
      cout<<"us["<<i<<"]: "<<us[i].transpose()<<endl;
      cout<<"ts["<<(i+1)<<"]: "<<ts[i+1]<<endl;
    }
  }
  else
  {
    double tf = 3;    // time horizon
    params.GetDouble("tf", tf);
    double h = tf/N;   // time step

    // initial state
    xs[0] = x0;
    us[0] = Vector2d(.5, .2);

    // Set times:
    for (int k = 0; k <=N; ++k)
      ts[k] = k*h;

    //////////////Creating the problem////////
    //Temporary point3d state:
    Point3dState projected_state;

    projectmanifold(xs[0],projected_state);

    gps(zs[0], ts[0], projected_state, us[0]);

    cout<<"Zs [0]: "<<zs[0].transpose()<<endl;

    //Using the controls Update the sys trajectory and find the sensor measurements:
    //Reset the system:
    sys.reset(xs[0], ts[0]);
    for(int i = 0; i< N; ++i)
    {
      //Set Control:
      if(i < N/2)
        us[i] = Vector2d(.5, .2);
      else
        us[i] = Vector2d(-.5, -.2);

      //Forward Step:
      sys.Step1(xs[i+1],us[i], h);
      //Project the state
      projectmanifold(xs[i+1],projected_state);

      cout<<"Xs["<<(i+1)<<"]: "<<xs[i+1].transpose()<<endl;
      //cout<<"Point3dState: "<<projected_state.q.transpose()<<endl;
      gps(zs[i+1], ts[i+1], projected_state, us[i]);
      //cout<<"Zs ["<<(i+1)<<"]: "<<zs[i+1].transpose()<<endl;
    }
  }
  RccarView view(sys, &xs);
  
  viewer->Add(view);
  
  //Assign  Zs:
  cost.SetReference(&zs, &mup);//Set reference for zs
  //Create Gauss newton estimation problem
  GnDoep<Vector4d, 4, 2, Dynamic, 9, Vector3d, 3, Point3dState, 6> gn(sys, gps, cost, ts, xs, us, p0, &projectmanifold);  
  getchar();


  struct timeval timer;
  /*
  cout<<"Varying Parameter and finding the cost and new state info: "<<endl;
  //Temporary sensor measurement:
  Vector3d ztemp;
  Vector8d residual;
  residual.resize(8);
  sys.reset(xs[0],ts[0]);
  //Testing Cost Function
  for(int i = 0; i< N; ++i)
  {
    //Forward Step:
    sys.Step3(xs[i+1],us[i], w0, h, &p0);
    //Project the state
    projectmanifold(xs[i+1],projected_state);

    cout<<"Xs["<<(i+1)<<"]: "<<xs[i+1].transpose()<<endl;
    //cout<<"Point3dState: "<<projected_state.q.transpose()<<endl;
    gps(ztemp, ts[i+1], projected_state, us[i]);
    //cout<<"Zs ["<<(i+1)<<"]: "<<ztemp.transpose()<<endl;
    //Cost Measurement:

    cout<<"Cost ["<<i<<"]: "<<cost.L(ts[i+1], ztemp, w0, p0, h)<< endl;
    cost.Res(residual, ts[i+1], ztemp, w0, p0, h);
    cout<<"Residual ["<<i<<"]: "<<residual.transpose()<< endl;

  }
  */

  gn.debug = false; // turn off debug for speed
  //getchar();

  for (int i = 0; i < iters; ++i) {
    timer_start(timer);
    gn.Iterate();
    long te = timer_us(timer);
    cout << "Parameter: "<< p0 << endl;
    cout << "Iteration #" << i << " took: " << te << " us." << endl;
    getchar();
  }


  //  xs[1][3]  velocity
  //atan(us[0][1]) steering angle
 
  cout << "done!" << endl;
  while(1)
    usleep(10);    
}


#define DISP

int main(int argc, char** argv)
{

  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("../../bin/gnrccarestimation.cfg"); 


#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "../../logs/rccar/frames/frame";

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
