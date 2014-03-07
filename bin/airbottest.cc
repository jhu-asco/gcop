#include <iomanip>
#include <iostream>
#include "dmoc.h"
#include "viewer.h"
#include "airbotview.h"
#include "utils.h"
#include "lqcost.h"
#include "params.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<MbsState> AirbotDmoc;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer) {
    viewer->SetCamera(45.25, 37, -0.14, 0.05, -1.75);
  }

  int N = 128;      // discrete trajectory segments
  double tf = 2;   // time-horizon

  params.GetInt("N", N);  
  params.GetDouble("tf", tf);  

  double h = tf/N; // time-step

  // system
  Airbot sys;


  VectorXd qv0(24);
  params.GetVectorXd("x0", qv0);  

  MbsState x0(7);
  SE3::Instance().rpyxyz2g(x0.gs[0], qv0.head(3), qv0.segment<3>(3));
    
  x0.r = qv0.segment<6>(6);
  x0.vs[0] = qv0.segment<6>(12);
  x0.dr = qv0.tail(6);
  sys.Rec(x0,h);


  VectorXd qvf(24);
  params.GetVectorXd("xf", qvf);  

  MbsState xf(7);
  SE3::Instance().rpyxyz2g(xf.gs[0], qvf.head(3), qvf.segment<3>(3));
    
  xf.r = qvf.segment<6>(6);
  xf.vs[0] = qvf.segment<6>(12);
  xf.dr = qvf.tail(6);


  LqCost<MbsState> cost(sys.X, (Rn<>&)sys.U, tf, xf);

  //  cost.Qf(3,3) = 50; cost.Qf(4,4) = 50; cost.Qf(5,5) = 50;
  //  cost.Qf(11,11) = 5; cost.Qf(12,12) = 5; cost.Qf(13,13) = 5;
  
  VectorXd Q(24);
  VectorXd R(10);  
  VectorXd Qf(24);

  params.GetVectorXd("Q", Q);
  params.GetVectorXd("R", R);
  params.GetVectorXd("Qf", Qf);

  cout << "Q=" << Q << endl;
 
  cost.Q = Q.asDiagonal();
  cost.R = R.asDiagonal();
  cost.Qf = Qf.asDiagonal();

  // times
  vector<double> ts(N+1);
  for (int k = 0; k <= N; ++k)
    ts[k] = k*h;

  // states
  vector<MbsState> xs(N+1, x0);

  double m = sys.links[0].m + sys.links[1].m  + sys.links[2].m + sys.links[3].m + sys.links[4].m  + sys.links[5].m + sys.links[6].m;
	cout<<"M = "<<m<<endl;
  
  // initial controls (e.g. hover at one place)
  VectorXd u(10);
  u.setZero();
  u(3) = 9.81*m;
  u(1) = .00;
  vector<VectorXd> us(N, u);
	cout<<"u[3]: "<<u[3]<<endl;

  AirbotView view(sys, &xs);
  if (viewer)
    viewer->Add(view);
  getchar();

  AirbotDmoc dmoc(sys, cost, ts, xs, us);
  params.GetDouble("mu", dmoc.mu);

  struct timeval timer;
  //  dmoc.debug = false; // turn off debug for speed

  Matrix4d g = Matrix4d::Identity();
  VectorXd gp(3);
  params.GetVectorXd("gp", gp);
  cout << "gp=" << gp << endl;
  g.topRightCorner<3,1>() = gp;
  CylinderView bv(.05, 1, &g);
  viewer->Add(bv);

	//Lets print stuff to compare:
	for(int count = 0;count<(sys.nb);count++)
	{
		cout<<"Ds["<<sys.links[count].name<<"]"<<endl<<sys.links[count].ds<<endl;
		cout<<"I["<<sys.links[count].name<<"]"<<endl<<sys.links[count].I<<endl;
		cout<<"M["<<sys.links[count].name<<"]"<<endl<<sys.links[count].m<<endl;

	}
	for(int count = 0;count<(sys.nb)-1;count++)
	{
		cout<<"Joint["<<sys.joints[count].name<<"].gc"<<endl<<sys.joints[count].gc<<endl;
		cout<<"Joint["<<sys.joints[count].name<<"].gp"<<endl<<sys.joints[count].gp<<endl;
		cout<<"Joint["<<sys.joints[count].name<<"].a"<<endl<<sys.joints[count].a<<endl;
	}
	
	cout<<"xf.r"<<endl<<xf.r<<endl;
	cout<<"Cost.Qf"<<endl<<cost.Qf<<endl;
	cout<<"Cost.Q"<<endl<<cost.Q<<endl;
	cout<<"Cost.R"<<endl<<cost.R<<endl;


  for (int i = 0; i < 50; ++i) {    
    timer_start(timer);
    dmoc.Iterate();
    long te = timer_us(timer);
    cout << "Iteration #" << i << ": took " << te << " us." << endl;        
        getchar();
  }

  int Nd;
  params.GetInt("Nd", Nd);
  vector<MbsState> xds(Nd+1, x0);  
  int d = N/Nd;
  for (int i=Nd, j = N; i>=0 && j>= 0; --i, j-=d)
    xds[i] = xs[j];

  //  xs = xds;
  
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
    params.Load("airflip.cfg");


#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "../../logs/airbot/frames/frame";

  pthread_t dummy;
  pthread_create( &dummy, NULL, (void *(*) (void *)) solver_process, viewer);

#else
  solver_process(0);
#endif


#ifdef DISP
  viewer->SetColor((float[3]){1,1,1});
  viewer->Start();
#endif


  return 0;
}
