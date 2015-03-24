#include <iomanip>
#include <iostream>
#include "viewer.h"
#include "rccarview.h"
#include "utils.h"
#include "rnlqcost.h"
#include "params.h"
#include "flatoutputtparam.h"

using namespace std;
using namespace Eigen;
using namespace gcop;


Params params;

void solver_process(Viewer* viewer)
{
  if (viewer)
    viewer->SetCamera(32, 47, -0.3, 0.4, -1.6);

  int N = 64;        // number of segments
  double tf = 5;    // time horizon

  int iters = 30;

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);

  params.GetInt("iters", iters);
  

  double h = tf/N;   // time step

  Rccar sys;

  //  sys.U.lb[1] = tan(-M_PI/5);
  //  sys.U.ub[1] = tan(M_PI/5);

  // initial state
  Vector4d x0(1,1,0,0);
  params.GetVector4d("x0", x0);

  // final state
  Vector4d xf(0,0,0,0);
  params.GetVector4d("xf", xf);  

  // cost
  RnLqCost<4, 2> cost(sys, tf, xf);
  VectorXd Q(4);
  if (params.GetVectorXd("Q", Q))
    cost.Q = Q.asDiagonal();
  
  VectorXd Qf(4);
  if (params.GetVectorXd("Qf", Qf))
    cost.Qf = Qf.asDiagonal();
  
  VectorXd R(2);
  if (params.GetVectorXd("R", R)) 
    cost.R = R.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <=N; ++k)
    ts[k] = k*h;

  // states
  vector<Vector4d> xs(N+1);
  // initial state
  xs[0] = x0;
	cout<<"N: "<<N<<endl;

  // initial controls
	//Update using initial controls:
  vector<Vector2d> us(N);
	sys.reset(xs[0],ts[0]);
  for (int i = 0; i < N; ++i) {
		if(i < N/2)
			us[i] = Vector2d(.01, .1);
		else
			us[i] = Vector2d(.01, -.1);
		sys.Step_internalinput(xs[i+1], us[i], ts[i+1]-ts[i]);
		cout<<"xs["<<i+1<<"]"<<xs[i+1].transpose()<<endl;//#DEBUG
		cout<<"us["<<i<<"]"<<us[i].transpose()<<endl;//#DEBUG
  }
  us[N-1] = Vector2d(0,0);
  
	RccarView view(sys, &xs);
  
  viewer->Add(view);
	getchar();

	//Create FlatOutputTparam<T,nx,nu>:
	int numberofknots = 10;
	FlatOutputTparam<Vector4d, 4, 2> fp(sys, 2, numberofknots, 2);//4 is the number of knots, 2 is the number of flat outputs, 2 is the number of derivatives needed by the system
	VectorXd s(numberofknots*2);//ny*numberofknots
	fp.To(s, ts, xs, us);
	getchar();
	//Using the parameters evaluate xs, us:
	fp.From(ts, xs,us, s);
	for(int count = 0; count < us.size(); count++)
	{
		cout<<"xs["<<count+1<<"]"<<xs[count+1].transpose()<<endl;//#DEBUG
		cout<<"us["<<count<<"]"<<us[count].transpose()<<endl;//#DEBUG
  }
	getchar();

  

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
    params.Load("../../bin/rccar.cfg"); 


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
