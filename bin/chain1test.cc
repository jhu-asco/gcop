#include <iomanip>
#include <iostream>
#include "utils.h"
#include "chain1.h"

using namespace std;
using namespace Eigen;
using namespace gcop;

//typedef Matrix<double, 8, 1> Vector8d;

void solver_process()
{

  int N = 200;     // discrete trajectory segments
  double tf = 1;   // time-horizon
  double h = (tf/(double)N);
	cout<<"h"<<h<<endl;


  // system
  Chain1 sys;
	//int nb = sys.nb;
	//print sys details:
	//Printing the mbsmodel params:
	for(int count = 0;count<(sys.nb);count++)
	{
		cout<<"Ds["<<sys.links[count].name<<"]"<<endl<<sys.links[count].ds<<endl;
		cout<<"I["<<sys.links[count].name<<"]"<<endl<<sys.links[count].I<<endl;
	}
	for(int count = 0;count<(sys.nb)-1;count++)
	{
		cout<<"Joint["<<sys.joints[count].name<<"].gc"<<endl<<sys.joints[count].gc<<endl;
		cout<<"Joint["<<sys.joints[count].name<<"].gp"<<endl<<sys.joints[count].gp<<endl;
		cout<<"Joint["<<sys.joints[count].name<<"].a"<<endl<<sys.joints[count].a<<endl;
	}

  //sys.debug = true;


  MbsState x(4);
  x.gs[0].setIdentity();
	x.vs[0].setZero();
	x.dr.setZero();
  x.r[0] = 3.14;
  x.r[1] = -1.57;
  x.r[2] = 1.57;
  sys.Rec(x,h);

  // states
  vector<MbsState> xs(N+1, x);

  // initial controls (e.g. hover at one place)
  VectorXd u(9);
  u.setZero();
	u[7] = 0.1;
  vector<VectorXd> us(N, u);

  struct timeval timer;
  //  dmoc.debug = false; // turn off debug for speed
	MatrixXd M(sys.nb + 5, sys.nb + 5);
	sys.Mass(M,xs[0]);
	cout<<"Mass matrix"<<endl<<M<<endl;

  for (int i = 0; i < N; ++i) {    
    timer_start(timer);
    sys.Step(xs[i+1], i*h, xs[i], us[i], h);
    long te = timer_us(timer);
		cout<<"xsnew.gs0"<<endl<<xs[i+1].gs[0]<<endl;
		cout<<"jointangles"<<xs[i+1].r<<endl;
		getchar();
    cout << "Iteration #" << i << ": took " << te << " us." << endl;
  }
  cout << "dr=" << xs[1].dr << endl;

  cout << "done!" << endl;

   while(1)
    usleep(10);    
}



int main(int argc, char** argv)
{

    solver_process();

  return 0;
}
