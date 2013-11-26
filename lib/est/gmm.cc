#include <assert.h>
#include <algorithm>
#include "utils.h"
#include "gmm.h"
//#include "itpp/stat/misc_stat.h"

using namespace gcop;
using namespace std;
//using namespace itpp;


Gmm::Gmm(int n, int k) :
  k(k), ns(k, Normal(n)), ws(k), cdf(k)
{
  assert(n > 0);
  assert(k > 0);

    //  ns = new Normal*[k];
    //  ws = new double[k];
    //  cdf = new double[k];

  for (int i = 0; i < k; ++i) {
    //    ns[i] = new Normal(n);
    ws[i] = 1/((double)k);
    cdf[i] = (i+1)/((double)k);

    ns[i].mu.setZero();
    ns[i].P = VectorXd::Constant(n, 1).asDiagonal();
  }

  tol = .01;
}

Gmm::~Gmm()
{
  /*
  for (int i = 0; i < k; ++i)
    delete ns[i];
  delete cdf;
  delete ws;
  delete ns;
  */
}


bool Gmm::Update()
{
  if (k == 1) {
    ws[0] = 1;
    cdf[0] = 1;
    return ns[0].Update();

  } else {
    /*
    bool ok = true;
    double wn = 0;
    for (int i = 0; i < k; ++i) {
      ok = ok && ns[i]->Update();
      wn += ws[i];
    }

    for (int i = 0; i < k; ++i) {
      ws[i] /= wn;
      cdf[i] = (i ? cdf[i - 1] + ws[i] : ws[i]); 
    }
    return ok;
    */    


    // No need to call update since it is called after EM

    bool ok = true;
    for (int i = 0; i < k; ++i) {
      cdf[i] = (i ? cdf[i - 1] + ws[i] : ws[i]); 
      ok = ok && ns[i].pd;
    }
    return ok;
  }
}


double Gmm::L(const VectorXd &x) const
{
  if (k == 1) {
    return ns[0].L(x);
  } else {
    
    double l = 0;
    for (int i = 0; i < k; ++i)
      l += ws[i]*ns[i].L(x);    
    return l;
  }
}


double Gmm::Sample(VectorXd &x)
{
  if (k == 1) {
    return ns[0].Sample(x);

  } else {
    // for now this is unefficient if k is big
    // TODO: implement as binary search
    double uc = rand()/(double)RAND_MAX;
    int i = 0;
    while (uc > cdf[i])
      ++i;
    
    assert(i < k);
    return ns[i].Sample(x);
  }
}

void Gmm::Init(const VectorXd &xlb, const VectorXd &xub)
{
  VectorXd dx = xub - xlb;   // range
  VectorXd r = dx/pow(k, 1.0/dx.size())/2;     // radius
  MatrixXd P = (r.cwiseProduct(r)).asDiagonal();
  for (int i = 0; i < k; ++i) {    
    ns[i].mu = xlb + dx.cwiseProduct(VectorXd::Random(xlb.size()));
    ns[i].P = P;
    ws[i] = 1.0/k;
    ns[i].Update();
  }

  Update();
}



void Gmm::Fit(const vector<pair<VectorXd, double> > &xps, int iter, const MatrixXd *S)
{
  int N = xps.size();

  if (k == 1) {
    ns[0].Fit(xps);
    if (S)
      ns[0].P += *S;
    return;
  }

  assert(N > 0);

  double ps[N][k];

  for (int l = 0; l < iter; ++l) {

    // E-step

    for (int j = 0; j < N; ++j) {
      
      const VectorXd &x = xps[j].first;      
      double p = xps[j].second;

      double norm = 0;
      double *psj = ps[j];

      for (int i = 0; i < k; ++i) {
        psj[i] = p*ns[i].L(x);    // likelihood of each sample        
        norm += psj[i];
      }

      //    assert(norm > 1e-10);
      //      cout << "    normalized: ps[" << j << "]=";
      for (int i = 0; i < k; ++i) {
        psj[i] /= norm;
        //        cout << psj[i] << " ";
      }    
      //      cout << endl;
    }  
    
    
    // M-step
    double maxd = 0;

    for (int i = 0; i < k; ++i) {
      double t1 = 0;
      VectorXd t2 = VectorXd::Zero(ns[0].mu.size());
      MatrixXd t3 = MatrixXd::Zero(ns[0].mu.size(), ns[0].mu.size());
      
      for (int j = 0; j < N; ++j) {
        const VectorXd &x = xps[j].first;
        t1 += ps[j][i];
        t2 += ps[j][i]*x;
        t3 += (ps[j][i]*x)*x.transpose();
      }
      
      ws[i] = t1;

      VectorXd mu = t2/t1;

      double d = (mu - ns[i].mu).norm();
      if (maxd < d)
        maxd = d;

      ns[i].mu = mu;
      ns[i].P = (t3 - t2*(t2.transpose()/t1))/t1;        
      if (S)
        ns[i].P += *S;


      if (!ns[i].Update()) // set Pinv, det, norm
        return;
    }
    
    if (maxd < tol) {
      //cout << "[W] Gmm::Fit: tolerance " << maxd << " reached after " << l << " iterations!" << endl;
      break;
    }
  }
}
