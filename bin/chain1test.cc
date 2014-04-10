#include <iomanip>
#include <iostream>
#include "dmoc.h"
#include "viewer.h"
#include "chainview.h"
#include "utils.h"
#include "lqcost.h"
#include "params.h"
#include "mbscontroller.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

typedef Dmoc<MbsState> ChainDmoc;

Params params;

void solver_process(Viewer* viewer)
{
  if (viewer) {
    viewer->SetCamera(45.25, 37, -0.14, 0.05, -1.75);
  }

  int N = 128;      // discrete trajectory segments
  double tf = 5;   // time-horizon
  int basetype = 1;
  const int FLOAT = 1;
  const int FIXED = 0;
  
  params.GetInt("N", N);
  params.GetDouble("tf", tf);  
  params.GetInt("fixed",basetype);
  assert((basetype == 0)||(basetype == 1));
  //basetype 1 is false that is floating
  //basetype 0 is true that is fixed
  

  int nb = 3;     // nof bodies
  params.GetInt("nb", nb);
  int n = nb -1;
  cout<<"n "<<n<<endl;
  
  double h = tf/N; // time-step
  
  // system
  Chain sys(nb,(basetype == FIXED));

  params.GetInt("method", sys.method);
  params.GetInt("iters", sys.iters);

  params.GetVectorXd("ulb", sys.U.lb);
  params.GetVectorXd("uub", sys.U.ub);


  // acceleration due to gravity
  params.GetVector3d("ag", sys.ag);
  VectorXd qv0(2*n+ 12*basetype);
  params.GetVectorXd("x0", qv0);  

  MbsState x0(nb,(basetype==0));
	if(basetype == FLOAT)
	{
		SE3::Instance().rpyxyz2g(x0.gs[0], qv0.head(3), qv0.segment<3>(3));
		x0.vs[0] = qv0.segment<6>(6*basetype+n);
		cout<<"1"<<endl;
	}
  x0.r = qv0.segment(6*basetype,n);
  x0.dr = qv0.tail(n);
  sys.Rec(x0, h);

  VectorXd qvf(2*n + 12*basetype);
  params.GetVectorXd("xf", qvf);  

  MbsState xf(nb,(basetype==0));
	if(basetype==FLOAT)
	{
		SE3::Instance().rpyxyz2g(xf.gs[0], qvf.head(3), qvf.segment<3>(3));
		xf.vs[0] = qvf.segment<6>(6*basetype+n);
		cout<<"2"<<endl;
	}
  xf.r = qvf.segment(6*basetype,n);
  xf.dr = qvf.tail(n);

  LqCost<MbsState> cost(sys.X, (Rn<>&)sys.U, tf, xf);
  
  VectorXd Q(2*n + 12*basetype);
  VectorXd R(n + 6*basetype);
  VectorXd Qf(2*n + 12*basetype);

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
  vector<MbsState> xs(N+1, x0);

  double m = sys.links[0].m + sys.links[1].m  + sys.links[2].m;
  cout<<"Mass: "<<m<<endl;
  
  // initial controls (e.g. hover at one place)
  VectorXd u(6*basetype + n);
  u.setZero();
	if(basetype==FLOAT)
		u(5) = -sys.ag[2]*m;
  vector<VectorXd> us(N, u);

  ChainView view(sys, &xs);
  if (viewer)
    viewer->Add(view);

  // @MK: this is the new part, initialize trajectory using a controller
  MbsController ctrl(sys);
  params.GetVectorXd("Kp", ctrl.Kp);
  params.GetVectorXd("Kd", ctrl.Kd);

  for (int i = 0; i < xs.size()-1; ++i) {
    double t = i*h;
    // ctrl.Set(us[i], t, xs[i]); 
    sys.Step(xs[i+1], i*h, xs[i], us[i], h);
  }

 cout << "dr=" << xs[1].dr << endl;
  // see the result before running optimization
  getchar();

  //ChainDmoc dmoc(sys, cost, ts, xs, us);
  //params.GetDouble("mu", dmoc.mu);

  //params.GetDouble("eps", dmoc.eps);

  //  dmoc.debug = false; // turn off debug for speed


  
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
    params.Load("../../bin/chainopt.cfg");


#ifdef DISP
  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "../../logs/chain/frames/frame";

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

