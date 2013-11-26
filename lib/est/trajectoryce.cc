#include "trajectoryce.h"
#include "utils.h"
#include <assert.h>

using namespace std;
using namespace dgc;
using namespace itpp;

TrajectoryCe::TrajectoryCe(const System &sys,
                           int sn,
                           int k,
                           int ni,
                           const int *is,
                           const double *inj) :
  Ce((ni > 0 ? ni : sys.n)*(sn-1), k, inj), si(sys), sf(sys), sn(sn), ni(ni), is(0), pad(false)
{
  assert(sn > 1);

  // set indices
  if (ni > 0) {
    // actual indices should be passed
    assert(is);
    this->is = new int[ni];
    memcpy(this->is, is, ni*sizeof(int));
  }
}


TrajectoryCe::~TrajectoryCe()
{
  delete[] is;
}


void TrajectoryCe::TrajToVec(double *z, const Trajectory &traj)
{
  for (int i = 1; i < sn; ++i) {
    if (ni > 0) {
      for (int j = 0; j < ni; ++j) {
        z[(i-1)*ni + j] = traj.states[MIN(i, traj.sn-1)]->x[is[j]];
      }
    } else {
      memcpy(z + (i-1)*traj.sys.n, traj.states[MIN(i, traj.sn-1)]->x, traj.sys.n*sizeof(double));
    }
  }
}

void TrajectoryCe::VecToTraj(Trajectory &traj, const double *z)
{
  assert(traj.sn >= sn);
  *traj.states[0] = si;
  *traj.states[sn] = sf;
  for (int i = 1; i < sn; ++i) {
    if (ni > 0)       
      for (int j = 0; j < ni; ++j)
        traj.states[i]->x[is[j]] = z[(i-1)*ni + j];
    else
      memcpy(traj.states[i]->x, z + (i-1)*traj.sys.n, traj.sys.n*sizeof(double));      
  }
  traj.SetTime(0, sf.t/sn);
}


void TrajectoryCe::AddSample(const Trajectory &traj, double c)
{
  //  assert(this->traj.sn <= traj.sn);

  vec Z(n);

  TrajToVec(Z._data(), traj);

  // when padding is on: sn >= traj.sn, otherwise sn <= traj.sn and MIN is unnecessary

  Ce::AddSample(Z, c);
}


double TrajectoryCe::Draw(Trajectory &traj)
{
  assert(traj.sn == sn);

  vec Z(n);
  double p = Ce::Draw(Z);
  
  VecToTraj(traj, Z._data());

  return p;
}


void TrajectoryCe::States(vector<State> &states, vector<double> &ps, int N, int sn)
{
  states.clear();
  ps.clear();

  if (sf.t < EPS)
    return;  

  Trajectory traj(si.sys, this->sn);  
  for (int j = 0; j < N; ++j) {
    double p = Draw(traj);   
    for (int i = 0; i <= sn; ++i) {
      State s(si.sys);
      s.t = si.t + RND*(sf.t - si.t);
      traj.Get(s);
      states.push_back(s);
      ps.push_back(p);
    }
  }
}
