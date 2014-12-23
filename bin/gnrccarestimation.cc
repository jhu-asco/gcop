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
  vector<Vector2d> us(N);
  //Sensor
  vector<Vector3d> zs(N/2);//Same as ts_sensor
  vector<double> ts_sensor(N/2);

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
 
  double tf = 3;    // time horizon
  params.GetDouble("tf", tf);
  double h = tf/N;   // time step

  // initial state
  xs[0] = x0;
  us[0] = Vector2d(.5, .2);

  // Set times:
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;
  //Set sensor times:
  for (int k = 0; k <(N/2); ++k)
    ts_sensor[k] = 2*k*h;

  //////////////Creating the problem////////
  //Temporary point3d state:
  Point3dState projected_state;

  //projectmanifold(xs[0],projected_state);

  //gps(zs[0], ts[0], projected_state, us[0]);

  //cout<<"Zs [0]: "<<zs[0].transpose()<<endl;

  //Using the controls Update the sys trajectory and find the sensor measurements:
  //Reset the system:
  sys.reset(xs[0], ts[0]);
  int sensor_index = 0;
  for(int i = 0; i< N; ++i)
  {
    //Set Control:
    if(i < N/2)
      us[i] = Vector2d(.5, .2);
    else
      us[i] = Vector2d(-.5, -.2);

    //Forward Step:
    sys.Step1(xs[i+1],us[i], h);

    cout<<"Xs["<<(i+1)<<"]: "<<xs[i+1].transpose()<<endl;
    if((ts_sensor[sensor_index] - ts[i])>= 0 && (ts_sensor[sensor_index] - ts[i+1]) < 0)
    {
      int near_index = (ts_sensor[sensor_index] - ts[i]) > -(ts_sensor[sensor_index] - ts[i+1])?(i+1):i;
      //Project the state
      projectmanifold(xs[near_index],projected_state);
      //cout<<"Point3dState: "<<projected_state.q.transpose()<<endl;
      gps(zs[sensor_index], ts[near_index], projected_state, us[near_index]);
      cout<<"Zs ["<<(sensor_index)<<"]: "<<zs[sensor_index].transpose()<<"ts_sensor: "<<ts_sensor[sensor_index]<<endl;
      sensor_index = sensor_index < (ts_sensor.size()-1)?sensor_index+1:sensor_index;
    }
  }
  RccarView view(sys, &xs);
  
  viewer->Add(view);
  getchar();
  
  //Assign  Zs:
  cost.SetReference(&zs, &mup);//Set reference for zs
  //Create Gauss newton estimation problem
  GnDoep<Vector4d, 4, 2, Dynamic, 9, Vector3d, 3, Point3dState, 6> gn(sys, gps, cost, ts, xs, us, p0, ts_sensor, &projectmanifold);  
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
