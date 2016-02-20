#include <iomanip>
#include <iostream>
#include "viewer.h"
#include "utils.h"
#include "qrotorsystemid.h"
#include "params.h"
#include <fstream>
#include <string>

using namespace std;
using namespace Eigen;
using namespace gcop;

Params params;

void solver_process()
{
   int N = 100;
   double tf = 2;

  QRotorIDState x0;
  vector<QRotorSystemIDMeasurement> systemid_measurements;
  //SO3
  SO3 &so3 = SO3::Instance();

  if(!params.Exists("dataFile"))
  {
    params.GetInt("N", N);
    params.GetDouble("tf", tf);
    double h = tf/N; // time-step
    // system
    QRotorIDModel sys;
    params.GetVector3d("kp",sys.kp);
    params.GetVector3d("kd",sys.kd);
    params.GetDouble("kt",sys.kt);
    params.GetVector3d("a0",sys.a0);
    params.GetVector3d("tau0",sys.tau0);

    x0.Clear();//Set to 0,0,0

    double thrust_bias = (9.81/sys.kt);
    double omega = 0.1;

    //Set SystemID Inputs
    systemid_measurements.resize(N);//control is 4x1 vector. First input is the thrust, second is the commanded rpydot

    QRotorIDState step, temp_step;//Iterating step
    step = x0;

    //Generate system trajectory:
    Vector3d current_rpy;
    for(int i = 0; i < N; i++)
    {
      double t = i*h;
      //Record current measurement and control,
      so3.g2q(current_rpy, step.R);
      systemid_measurements[i].position = step.p;
      systemid_measurements[i].rpy = current_rpy;
      systemid_measurements[i].t = t;
      Vector4d &current_control = systemid_measurements[i].control;
      current_control<<(thrust_bias + 1*cos(omega*t)), 1*cos(omega*t), 1*sin(omega*t), 1*sin(omega*t);
      cout<<i<<" "<<step.p.transpose()<<" "<<current_rpy.transpose()<<" "<<step.u.transpose()<<" "<<current_control.transpose()<<endl;
      //Propagate  the step
      sys.Step(temp_step,t,step,current_control,h,0,0,0,0);
      step = temp_step;
    }
    x0.p = x0.p;// + Vector3d(0.01, 0.01, 0.01);//Adding random noise to initial state for optimization;
  }
  else
  {
      string dataFileName;
      params.GetString("dataFile",dataFileName);
      cout<<"Find Data File: "<<dataFileName<<endl;
      //Open Data File:
      ifstream ifile(dataFileName);
      if(!ifile.is_open())
      {
          cerr<<"File Not Open"<<endl;
          return;
      }
      std::string temp;
      while(!ifile.eof())
      {
        //Read data into systemid_measurements:
        QRotorSystemIDMeasurement measurement;
        if(!getline(ifile,temp))//Get string until \n
          break;

        stringstream ss(temp);//Make string stream out of the temp string
        ss>>measurement.t;
        for(int i = 0; i < 3; i++)
        {
          ss>>measurement.position[i];
        }
        for(int i = 0; i < 3; i++)
        {
          ss>>measurement.rpy[i];
        }
        for(int i = 0; i < 4; i++)
        {
          ss>>measurement.control[i];
        }
        
        systemid_measurements.push_back(measurement);
        cout<<measurement.t<<" "<<measurement.position.transpose()<<" "<<measurement.rpy.transpose()<<" "<<measurement.control.transpose()<<endl;
      }
      //Get x0:
      x0.Clear();
      x0.p = systemid_measurements[0].position;
      const Vector3d &rpy = systemid_measurements[0].rpy;
      so3.q2g(x0.R,rpy);
      x0.u<<0,0,rpy(2);//0,0,yaw
  }

  Matrix7d stdev_gains;
  Vector6d mean_offsets;
  Matrix6d stdev_offsets;

  //Create QRotorSystemID:
  QRotorSystemID system_id;
  params.GetDouble("offsets_period",system_id.offsets_timeperiod);
  system_id.EstimateParameters(systemid_measurements, x0, &stdev_gains, &mean_offsets, &stdev_offsets);
}


int main(int argc, char** argv)
{

  if (argc > 1)
    params.Load(argv[1]);
  else
    params.Load("../../bin/qrotorsystemidtest.cfg");

  solver_process();


  return 0;
}
