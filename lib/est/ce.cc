#include <assert.h>
#include <algorithm>
#include "utils.h"
#include "ce.h"
#include <limits>

using namespace std;
using namespace gcop;


Ce::Ce(int n,
       int k,
       const MatrixXd *S) :
  n(n),
  gmm(n, k),
  S(S ? *S : MatrixXd::Zero(n,n)),
  rho(.1),
  alpha(.9),
  mras(false),
  b(1)
{
}


Ce::~Ce()
{ 
}


void Ce::Reset()
{
  zps.clear();
  cs.clear();
}

void Ce::AddSample(const VectorXd &z, double c) 
{ 
  zps.push_back(make_pair(z, c));
  cs.push_back(c);
}

bool zpSort(const pair<VectorXd, double> &zpa, const pair<VectorXd, double> &zpb)
{
  return zpa.second < zpb.second;
}

void Ce::Select()
{
  if (mras)
    return;
  
  int N = zps.size();
  int Ne = MIN((int)ceil(N*rho), N);
  
  if (Ne < N) {
    std::sort(zps.begin(), zps.end(), zpSort);
    std::sort(cs.begin(), cs.end());   // this is redundant but is kept for consistency
    zps.resize(Ne);
    cs.resize(Ne);
  }
}


bool Ce::Fit()
{       
  int N = zps.size();
  
  if (mras) {
    double cn = 0;
    for (int j = 0; j < N; ++j) {
      pair<VectorXd, double> &zp = zps[j];
      zp.second = exp(-b*cs[j]);     // cost
      cn += zp.second;               // normalizer
    }
    assert(cn > 0);
    for (int j = 0; j < N; ++j) 
      zps[j].second /= cn;               // normalize
    
    gmm.Fit(zps, alpha, 50, &S);
  } else {
    for (int j = 0; j < N; ++j) 
      zps[j].second = 1.0/N;             // probability
    
    gmm.Fit(zps, alpha, 50, &S);
  }

  return gmm.Update();
}

double Ce::Sample(VectorXd &z) 
{
  return gmm.Sample(z);
}

const VectorXd& Ce::Best() 
{
  return zps[0].first;
}
