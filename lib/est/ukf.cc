#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include "ukf.h"
#include "quat.h"

using namespace gcop;
using namespace std;
using namespace Eigen;

UKF::UKF(Model &model) :
  Filter(model), 
  Xs(new VectorXd*[2*model.nr + 1]),
  Xps(new VectorXd*[2*model.nr + 1]),
  Zs(new VectorXd*[2*model.nr + 1]),
  Ws(2*model.nr + 1),
  Wc(2*model.nr + 1),
  Pzz(model.nz, model.nz),
  Pxz(model.nr, model.nz),
  A(model.nr, model.nr)
{
  this->L = model.nr;
  this->a = .001;
  this->k = 0;
  this->b = 2;
  this->l = a*a*(L+k)-L;
  
  for (int i = 0; i < 2*L + 1; ++i) {    
    
    Xs[i] = new VectorXd(model.nx);
    Xps[i] = new VectorXd(model.nx);
    Zs[i] = new VectorXd(model.nz);
    
    if (i == 0) {
      Ws[0] = l/(L+l);
      Wc[0] = l/(L+l) + (1-a*a+b);
    } else {
      Ws[i] = 1/(2*(L+l));
      Wc[i] = 1/(2*(L+l));
    }
  }
}


UKF::~UKF()
{
  for (int i = 0; i < 2*L + 1; ++i) {
    delete Xps[i];
    delete Xs[i];
    delete Zs[i];
  }
  delete[] Xps;
  delete[] Xs;
  delete[] Zs;
}


bool UKF::Predict(const VectorXd *u,
                  bool cov)
{ 
  A = this->P.llt().matrixL();
  bool pd = true;

  if (!pd) {
    cout << "[W] UKF::Predict: Cholesky failed!" << endl;
    return false;
  }
  
  Points(Xs, x, sqrt(L+l)*A.transpose());

  x.setZero();  // zero mean
  
  for (int i = 0; i < 2*L + 1; ++i) {
    
    if (u)
      model.f(*Xps[i], *Xs[i], *u);
    else
      model.f(*Xps[i], *Xs[i]);

    x = x + Ws[i]*(*Xps[i]);
  }

  P.setZero(); // zero covariance

  if (model.quat) {
    
    double xn = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
    x[0] /= xn;
    x[1] /= xn;
    x[2] /= xn;
    x[3] /= xn;
    
    // inverse
    Quat qi(x.data());
    qi.Invert();

    VectorXd dx(L);
    for (int i = 0; i < 2*L+1; ++i) {
      // difference on the group
      VectorXd &xc = *Xps[i];
      Quat qc(xc.data());
      Quat dq = qi*qc;
      dq.Normalize();
      dq.ToExp(dx.data());
      dx.tail(dx.size() - 3) = xc.tail(xc.size() - 4) - x.tail(x.size() - 4);
      P = P + Wc[i]*dx*dx.transpose();
    }
  } else {
    for (int i = 0; i < 2*L+1; ++i) {
      VectorXd dx = (*Xps[i]) - x;
      P = P + Wc[i]*dx*dx.transpose();
    }
  }

  P = P + model.Q;
  return true;
}



bool UKF::Update(const VectorXd &z, bool cov)
{
  A = this->P.llt().matrixL();

  bool pd = true;
  if (!pd) {
    cout << "[W] UKF::Update: Cholesky failed!" << endl;
    return false;
  }

  Points(Xs, x, sqrt(L+l)*A.transpose());

  VectorXd xm = VectorXd::Zero(x.size());
  VectorXd zm = VectorXd::Zero(z.size());
  VectorXd zt(z.size());  // temp measurement

  for (int i = 0; i < 2*L+1; ++i) {
    VectorXd& xc = *Xs[i];
    VectorXd& zt = *Zs[i];
    model.h(zt, xc);
    zm = zm + Ws[i]*zt;
    xm = xm + Ws[i]*xc;
  }

  // normalize it
  if (model.quat) {
    zm.head(3).normalize();
    xm.head(4).normalize();
  }
  
  VectorXd dx(L);

  Pzz.setZero();
  Pxz.setZero();

  for (int i = 0; i < 2*L+1; ++i) {
    VectorXd dz = (*Zs[i]) - zm;
    if (model.quat) {
      VectorXd& xc = *Xs[i];
      Quat qmi(xm.data());
      Quat q(xc.data());
      qmi.Invert();
      Quat dq = qmi*q;
      dq.Normalize();
      dq.ToExp(dx.data());
      dx.tail(dx.size()-3) = xc.tail(xc.size()-4) - xm.tail(xm.size()-4); 
    } else {
      dx = (*Xs[i]) - xm;
    }
    Pzz = Pzz + Wc[i]*dz*dz.transpose();
    Pxz = Pxz + Wc[i]*dx*dz.transpose();
  }
  
  Pzz = Pzz + model.R;
  MatrixXd K = Pxz*Pzz.inverse();
  
  dx = K*(z - zm);

  if (model.quat) {
    Quat q(x.data());
    Quat dq;
    dq.FromExp(dx.data());
    q = q*dq;    
    q.Q(x.data());
    x.tail(x.size()-4) = x.tail(x.size()-4) + dx.tail(dx.size() - 3);
  } else {
    x = x + dx;
  }
  
  P = P - K*Pzz*K.transpose();
  
  return true;
}


void UKF::Points(VectorXd **Xs, 
                 const VectorXd &x,
                 const MatrixXd &A)
{
  (*Xs[0]) = x;


  for (int i = 0; i < L; ++i) {
    const VectorXd& dx = A.col(i);

    if (model.quat) {
      Quat q(x.data());
      Quat qw;

      VectorXd& xcp = *Xs[i + 1];
      qw.FromExp(dx.data());
      qw = q*qw;
      qw.Q(xcp.data());
      xcp.tail(xcp.size()-4) = x.tail(x.size()-4) + dx.tail(dx.size()-3);

      VectorXd& xcm = *Xs[i + 1 + L];
      Vector3d wm = -dx.head(3);
      qw.FromExp(wm.data());
      qw = q*qw;
      qw.Q(xcm.data());
      xcm.tail(xcm.size()-4) = x.tail(x.size()-4) - dx.tail(dx.size()-3);

    } else {
      (*Xs[i + 1]) = x + dx;
      (*Xs[i + 1 + L]) = x - dx;
    }
  }
}
