//#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <iostream>
//#include <Eigen/Sparse>
//#include <unsupported/Eigen/src/SparseExtra/UmfPackSupport.h>
#include "gsl/gsl_cdf.h"
#include "utils.h"

#include "gp.h"

using namespace gcop;
using namespace std;

using namespace Eigen;

//#define SYM


int main(int argc, char** argv)
{
  GP gp(2, 0);
  gp.sigma = .1;

  //  gp.Sample();
  
  //  gp.Train();
  
  for (int i = 0; i < 500; ++i) {
    VectorXd x(2);
    x(0) = 5*(RND - .5);
    x(1) = 5*(RND - .5);
    gp.Add(x, x.norm() - 1);
  }

  gp.OptParams();
  cout << "l=" << gp.l << endl;  

  int n = 10;
  MatrixXd ms(n,n);
  MatrixXd ss(n,n);

  double s = 0;

  VectorXd x(2);
  for(int i = 0; i < n; ++i) {
    x(0) = -2 + ((double)i)/(n-1)*4;
    for(int j = 0; j < n; ++j) {
      x(1) = -2 + ((double)j)/(n-1)*4;
      
      ms(j,i) = gp.Predict(x, &s);
      ss(j,i) = s;
    }
  }

  //  cout << gp.Xs << endl;
  //  cout << gp.fs << endl;

  cout << ms << endl;
  cout << endl;
  cout << ss << endl;

  return 1;

}
