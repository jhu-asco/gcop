#include <iostream>
#include "utils.h"
#include "gp.h"

using namespace gcop;
using namespace Eigen;

GP::GP(int d, int n) : 
  d(d), n(n), 
  Xs(d,n), fs(n), K(n,n), Ki(n,n), L(n,n), a(n),
  l(1), s(1), sigma(0), cf(true), eps(-1)
{
}

GP::~GP() 
{
}


void GP::Sample()
{
  for (int j = 0; j < Xs.cols(); ++j) {
    for (int i = 0; i < 2; ++i) {
      Xs(i,j) = 5*(RND - .5);
    }
    fs[j] = Xs.col(j).norm() - 1;
  }
}



bool GP::Add(const VectorXd &x, double f)
{
  int n = Xs.cols();
  
  // current data correlation vector
  VectorXd k(n);
  
  for (int i = 0; i < n; ++i)  {
    if (eps > 0 && (x - Xs.col(i)).norm() < eps)
      return false;
    k(i) = SqExp(Xs.col(i), x);
  }
  
  double kn = SqExp(x, x);
  if (sigma > 0)
    kn += sigma*sigma;
      
  // update covariance matrix
  K.conservativeResize(n + 1, n + 1);
  if (n>0) {
    K.block(0, n, n, 1) = k;
    K.block(n, 0, 1, n) = k.transpose();
  }
  K(n, n) = kn;
  
  // add values to value vector
  fs.conservativeResize(n + 1);
  fs[n] = f;

  // add data to data matrix
  Xs.conservativeResize(d, n + 1);
  Xs.col(n) = x;
  
  if (cf) {
    // update L
    VectorXd l = L.triangularView<Eigen::Lower>().solve(k);
    
    L.conservativeResize(n + 1, n + 1);
    if (n > 0) {
      L.block(0, n, n, 1).setZero();
      L.block(n, 0, 1, n) = l.transpose();
    }
    L(n, n) = sqrt(kn - l.dot(l));
    
    a = L.transpose().triangularView<Eigen::Upper>().solve(L.triangularView<Eigen::Lower>().solve(fs));

    n++;

  } else {
    // update K^{-1}
 
    double mu = 1/(kn - k.dot(Ki*k));
    VectorXd em = -mu*(Ki*k);
    MatrixXd M = Ki - (mu*k)*k.transpose();
    
    Ki.resize(n + 1, n + 1);
    Ki.block(0, 0, n, n) = M;
    Ki.block(0, n, n, 1) = em;
    Ki.block(n, 0, 1, n) = em.transpose();
    Ki(n, n) = mu;
    
    a = Ki.transpose()*fs;
  }

  return true;
}


void GP::Train(const MatrixXd &Xs, const VectorXd &fs)
{
  this->Xs = Xs;
  this->fs = fs;
  K.resize(Xs.cols(), Xs.rows());
  Ki.resize(Xs.cols(), Xs.rows());
  Train();
}


void GP::Train()
{
  //  MatrixXd K(n,n);
  n = Xs.cols();
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) { 
      if (i==j && sigma > 0)
        K(i,j) = SqExp(Xs.col(i), Xs.col(j)) + sigma*sigma;
      else
        K(i,j) = SqExp(Xs.col(i), Xs.col(j));
    }
  }

  LLT<MatrixXd> lltOfA(K); // compute the Cholesky decomposition of A
  L = lltOfA.matrixL();

  //  a = L.triangularView<Eigen::Lower>().solveInPlace(f);

  a = lltOfA.solve(fs);

  if (!cf)
    Ki = K.inverse();

  //  MatrixXd L = trans(chol(K));
  //  a = trans(S.L)\(S.L\S.fs);


  // S.lp = -S.fs'*S.a - sum(log(diag(S.L))) - S.n/2*log(2*pi);
  
  //  a = inv(sympd(K))*fs;


  //\  lp = -dot(fs, a) - sum(log(diag(L))) - n/2.0*log(2*M_PI);
 
}


double GP::LogL(double dll[2])
{
  if (dll) {
    MatrixXd dKdl(n,n); // der wrt l
    MatrixXd dKds(n,n); // der wrt s
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        VectorXd d = Xs.col(i) - Xs.col(j);
        double dn = d.dot(d);
        double se = s*s*exp(-dn/(2*l*l));
        dKdl(i,j) = se*dn/(l*l*l);
        dKds(i,j) = se*2/s;
      }
    }

    dll[0] = (a*a.transpose()*dKdl - L.transpose().triangularView<Eigen::Upper>().solve(L.triangularView<Eigen::Lower>().solve(dKdl))).trace()/2;
    
    dll[1] = (a*a.transpose()*dKds - L.transpose().triangularView<Eigen::Upper>().solve(L.triangularView<Eigen::Lower>().solve(dKds))).trace()/2;    
  }
  
  double w = 0;
  for (int i = 0; i < n; ++i)
    w += log(L(i,i));
  
  return -w - fs.dot(a) - n*log(2*M_PI)/2;
}


double GP::OptParams()
{
  double llmax = 0;
  double lmax = l;
  double dl = .01;
  for (l = dl; l < 5; l += dl) {
    Train();
    double ll = LogL();
    if (llmax < ll) {
      llmax = ll;
      lmax = l;
    }
  }
  l = lmax;
  return l;
}


double GP::Predict(const VectorXd &x, double *s) const
{
  double m = 0;
  
  VectorXd k(Xs.cols());

  for (int j = 0; j < Xs.cols(); ++j) {
    k(j) = SqExp(Xs.col(j), x);
    m += a(j)*k(j);
  }
  
  if (s) {
    if (cf) {
      VectorXd v = L.triangularView<Eigen::Lower>().solve(k);
      *s = SqExp(x, x) - v.dot(v);
    } else {
      *s = SqExp(x, x) - k.dot(Ki*k);
    }
    
    if (sigma > 0)
      *s += sigma*sigma;
  }

  return m;
}



double GP::SqExp(const VectorXd &xa, const VectorXd &xb) const
{
  VectorXd d = xa - xb;
  return s*s*exp(-d.dot(d)/(2*l*l));
}


double GP::PI(const VectorXd &x, double fmin) const
{
  double s;
  double m = Predict(x, &s);
  return ncdf((fmin - m)/s);
}
