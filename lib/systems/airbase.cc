#include "airbase.h"

using namespace gcop;
using namespace Eigen;

Airbase::Airbase(int nb,int j) : 
  Mbs(nb, 4+j) {
}

void Airbase::Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
                   MatrixXd *A, MatrixXd *B) 
{
  assert(f.size() > 6);
  assert(u.size() > 4);
  assert(u.size()-4 == f.size() -6);

  f[0] = u[0];
  f[1] = u[1];
  f[2] = u[2];
  f[3] = 0;
  f[4] = 0;
  f[5] = u[3];
  f.tail(f.size()-6) = u.tail(u.size()- 4);

  /*
  f[6] = u[4];
  f[7] = u[5];
  f[8] = u[6];
  f[9] = u[7];
  f[10] = u[8];
  f[11] = u[9];
  */

  /*f[6] = u[4] - .1*x.dr[0];
  f[7] = u[5] - .1*x.dr[1];
  f[8] = u[6] - .1*x.dr[2];
  f[9] = u[7] - .1*x.dr[3];
  f[10] = u[8] - .1*x.dr[4];
  f[11] = u[9] - .1*x.dr[5];
	*/

  if (A)
    A->setZero();

  if (B) {
    B->setZero();
    (*B)(0,0) = 1;
    (*B)(1,1) = 1;
    (*B)(2,2) = 1;

    for(int count = 5;count < f.size();count++)
      (*B)(count,count-2) = 1;
    cout<<"F.size: "<<f.size()<<endl;

    /*(*B)(7,5) = 1;
      (*B)(8,6) = 1;
    (*B)(9,7) = 1;
    (*B)(10,8) = 1;
    (*B)(11,9) = 1;
		*/
  }
}
