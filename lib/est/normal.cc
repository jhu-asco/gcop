#include <assert.h>
#include <algorithm>
#include "utils.h"
#include "normal.h"

using namespace gcop;
using namespace std;

Normal::Normal(int n) : 
  mu(VectorXd::Zero(n)),
  P(MatrixXd::Zero(n,n)),
  det(0),
  Pinv(MatrixXd::Zero(n,n)),
  pd(false),
  A(MatrixXd::Zero(n,n)),
  rn(VectorXd::Zero(n)),
  norm(0),
  bd(0)
{
}


Normal::Normal(const VectorXd &mu, const MatrixXd &P):
  mu(mu),
  P(P),
  det(0),
  Pinv(MatrixXd::Zero(mu.size(),mu.size())),
  pd(false),
  A(MatrixXd::Zero(mu.size(),mu.size())),
  rn(VectorXd::Zero(mu.size())),
  norm(0),
  bd(0)
{
  Update();
}


Normal::~Normal()
{
}


double Normal::L(const VectorXd &x) const
{
  if (!pd) {
    cout << "[W] Normal::L: not positive definite!" << endl;
  }
  
  VectorXd d = x - mu;
  return exp(-d.dot(Pinv*d))/2/norm;
}

bool Normal::Update()
{
  llt.compute(P);
  
  if (llt.info() == Eigen::Success) {
    A = llt.matrixL();
    Pinv = P.inverse();
    det = P.determinant();
    norm = sqrt(pow(2*M_PI, mu.size())*det);    
    pd = true;
  } else {
    cout << "[W] Normal::Update: cholesky failed: P=" << P << endl;
    pd = false;    
  }

  return pd;
}


double Normal::Sample(VectorXd &x)
{
  double p = 1;
  for (int i = 0; i < rn.size(); ++i) {
    rn(i) = random_normal();
    p *= rn(i);
  }
  
  x = mu + A*rn;
  return p;
}


void Normal::Fit(const vector<pair<VectorXd, double> > xws)
{
  int N = xws.size();

  mu.setZero();
  for (int j = 0; j < N; ++j) {
    const pair<VectorXd, double> &xw = xws[j];
    mu += xw.first*xw.second;
  }
  
  P.setZero();
  for (int j = 0; j < N; ++j) {
    const pair<VectorXd, double> &xw = xws[j];
    VectorXd dx = xw.first - mu;
    if (!bd) {
      P += (xw.second*dx)*dx.transpose();
    } else {
      int b = dx.size()/bd;
      for (int i = 0; i < b; ++i) {
        int bi = i*bd;
        VectorXd bdx = dx.segment(bi, bd);        
        P.block(bi, bi, bd, bd) += (xw.second*bdx)*bdx.transpose();       
      }
    }
  }
}
