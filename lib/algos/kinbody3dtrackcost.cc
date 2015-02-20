#include <assert.h>
#include "kinbody3dtrackcost.h"
#include "kinbody3dmanifold.h"
#include "rn.h"
#include "so3.h"

using namespace std;
using namespace gcop;
using namespace Eigen;

Kinbody3dTrackCost::Kinbody3dTrackCost(double tf, const Kinbody3dTrack &pg) :
  Cost(pg.sys, tf), pg(pg), useFdHess(false)
{
  
}


static Matrix3d r3hat(const Vector3d &a)
{
  SO3& so3 = SO3::Instance();
  Matrix3d a_hat;

  so3.hat(a_hat, a);
  //a_hat << 0,  -a(2),  a(1),
  //       a(2),    0,  -a(0),
  //      -a(1),  a(0),    0;
  return a_hat;
}


static Vector3d cross3(Vector3d a, Vector3d b)
{
  return a.cross(b);
}

void Kinbody3dTrackCost::pokeX(Matrix4d &xb, const Matrix4d &xa, const int i, const double eps)
{
  Kinbody3dManifold& km = Kinbody3dManifold::Instance();
  xb.setZero();
  Vector6d vh;
  vh.setZero();
  vh(i) = eps;
  km.Retract(xb, xa, vh);
} 

void Kinbody3dTrackCost::pokeX2(Matrix4d &xb, const Matrix4d &xa, const int i, const int j, const double epsi, const double epsj)
{
  Kinbody3dManifold& km = Kinbody3dManifold::Instance();
  xb.setZero();
  Vector6d vh;
  vh.setZero();
  vh(i) = epsi;
  vh(j) = epsj;
  km.Retract(xb, xa, vh);
}


Matrix6d Kinbody3dTrackCost::fdLxx(double t, double h, const Matrix4d &x,  const VectorXd *p)
{

  Matrix6d Lxx;
  Lxx.setZero();

  //cout << x << endl;
  //cout << *p << endl;

  //cout.precision(15);
  const double eps = 1e-7;
  for(int i = 0; i < 6; i++)
  {
    for(int j = 0; j < 6; j++)
    {
      if(i == j)
      {
         Matrix4d xp1, xp2, xp3, xp4, xp5;
         pokeX(xp1, x, j, 2*eps);
         pokeX(xp2, x, j, eps); 
         xp3 = x;
         pokeX(xp4, x, j, -eps);
         pokeX(xp5, x, j, -2*eps);

         double c1, c2, c3, c4, c5;
         c1 = obsCost(t, h, xp1, p);
         c2 = obsCost(t, h, xp2, p);
         c3 = obsCost(t, h, xp3, p);
         c4 = obsCost(t, h, xp4, p);
         c5 = obsCost(t, h, xp5, p);
         //cout << i << " " << j << endl;
         //cout << c1 << endl;
         //cout << c2 << endl;
         //cout << c3 << endl;
         //cout << c4 << endl;
         //cout << c5 << endl;
         //cout << (-c1 + 16.*c2 - 30.*c3 + 16.*c4 - c5)/(12.*eps*eps) << endl;
         Lxx(i,j) = (-c1 + 16.*c2 - 30.*c3 + 16.*c4 - c5)/(12.*eps*eps);
      }
      else
      {
         Matrix4d xp1, xp2, xp3, xp4;
         pokeX2(xp1, x, i, j, eps, eps);
         pokeX2(xp2, x, i, j, eps, -eps); 
         pokeX2(xp3, x, i, j, -eps, eps); 
         pokeX2(xp4, x, i, j, -eps, -eps); 

         double c1, c2, c3, c4, c5;
         c1 = obsCost(t, h, xp1, p);
         c2 = obsCost(t, h, xp2, p);
         c3 = obsCost(t, h, xp3, p);
         c4 = obsCost(t, h, xp4, p);
         Lxx(i,j) = (c1 - c2 - c3 + c4)/(4.*eps*eps);
      }
    }
  }
  return Lxx;
}

double Kinbody3dTrackCost::obsCost(double t, double h, const Matrix4d &x, const VectorXd *p)
{
  const Matrix3d &R = x.block<3,3>(0,0); // orientation
  const Vector3d &xp = x.block<3,1>(0,3);     // position
  int N = pg.Is.size() - 1;
  double L = 0;
  
  h = this->tf/N;

  int k = (int)round(t/h);
  assert(k >=0 && k <= N);
  const vector< pair<int, Vector3d> > &I = pg.Is[k];

  int i0 = 3*pg.extforce;
  for (int i = 0; i < I.size(); ++i) {
    int l = pg.cis[I[i].first];  // feature index
    const Vector3d &z = I[i].second;

    const Vector3d &pf = p->segment<3>(i0 + 3*l);
    Vector3d y = R.transpose()*(pf - xp);
    Vector3d r = (y - z)/pg.cp;

    L += r.dot(y - z)/2;         // feature cost
  }
  return L;
}

double Kinbody3dTrackCost::L(double t, const Matrix4d &x, const Vector6d &u, 
                          double h,
                          const VectorXd *p,
                          Vector6d *Lx, Matrix6d *Lxx,
                          Vector6d *Lu, Matrix6d *Luu,
                          Matrix<double, 6, 6> *Lxu,
                          VectorXd *Lp, MatrixXd *Lpp, MatrixX6d *Lpx)
{
  double L = 0;

  //cout << "x: " << endl << x << endl;

  if (Lu)
    Lu->setZero();

  if (Luu)
    Luu->setZero();

  if (Lxu)
    Lxu->setZero();
  
  if (Lx)
    Lx->setZero();

  if (Lxx)
    Lxx->setZero();

  if (Lp)
    Lp->setZero();

  if (Lpp)
    Lpp->setZero();

  if (Lpx)
    Lpx->setZero();
  
  //const Matrix3d &g = x.first;

  const Matrix3d &R = x.block<3,3>(0,0); // orientation
  const Vector3d &xp = x.block<3,1>(0,3);     // position
  
  int N = pg.Is.size() - 1;
  
  h = this->tf/N;

  //  int k = (h < 1e-16 ? N : (int)round(t/h));

  int k = (int)round(t/h);
  assert(k >=0 && k <= N);
  const vector< pair<int, Vector3d> > &I = pg.Is[k];

  //cout << "Kinbody3dtrackCost: k=" << k << " " << I.size() << endl;

  int i0 = 3*pg.extforce;
  for (int i = 0; i < I.size(); ++i) {
    int l = pg.cis[I[i].first];  // feature index
    const Vector3d &z = I[i].second;

    const Vector3d &pf = p->segment<3>(i0 + 3*l);
    Vector3d y = R.transpose()*(pf - xp);
    Vector3d r = (y - z)/pg.cp;

    L += r.dot(y - z)/2;         // feature cost

    //cout << "feature cost " << i << "=" << r.dot(y - z)/2 << endl;

    if (Lx) {
      //Lx->segment<3>(0) = Lx->segment<3>(0) - cross3(y, r);
      Lx->segment<3>(0) = Lx->segment<3>(0) - r3hat(y)*r;
      Lx->segment<3>(3) = Lx->segment<3>(3) - r;
    }

    if (Lxx) {
      if(!useFdHess)
      {
        Lxx->block<3,3>(0,0) += (r3hat(y).transpose()*r3hat(y))/pg.cp + r3hat(r)*r3hat(y)/2. + (r3hat(r)*r3hat(y)).transpose()/2.;
        Lxx->block<3,3>(0,3) += (-r3hat(r)/2. + r3hat(y)/pg.cp);
        Lxx->block<3,3>(3,0) += (-r3hat(r)/2. + r3hat(y)/pg.cp).transpose();      
        Lxx->block<3,3>(3,3) += Matrix3d::Identity()/pg.cp;
      }
      //*Lxx = Eigen::MatrixXd::Identity(6,6);
    }

    if (Lp)
      Lp->segment<3>(i0 + 3*l) = R*r;

    if (Lpp)
    {
      Lpp->block<3,3>(i0 + 3*l, i0 + 3*l) = Matrix3d::Identity()/pg.cp;
      //Lpp->block<3,3>(i0 + 3*l, i0 + 3*l) = Matrix3d::Identity();
    }

    if (Lpx) {
      Lpx->block<3,3>(i0 + 3*l, 0) = R*((r3hat(y)/pg.cp) - r3hat(r));
      Lpx->block<3,3>(i0 + 3*l, 3) = -R/pg.cp;
      //Lpx->block<3,3>(i0 + 3*l, 3) = Eigen::MatrixXd::Zero(3,3);
    }
  }

  if(Lxx)
  {
    if(useFdHess)
    {
      *Lxx = fdLxx(t, h, x, p);
    }
  }
  //cout << "t=" << t << " feature L=" << L << endl;

  if (pg.forces && k < N) {
    Vector6d du = u - pg.uos[k];
    Vector6d wdu = du.cwiseQuotient(pg.cw);

    L += du.dot(wdu)/2;
    if (Lu)
      *Lu = wdu;
    if (Luu)
      *Luu = Vector6d::Ones().cwiseQuotient(pg.cw).asDiagonal();
  }

  /*
  //if(t < this->tf - 1e-10)
  {
    if(Lu)
      *Lu *= h;
    if(Luu)
      *Luu *= h;
    if(Lx)
      *Lx *= h;
    if(Lxx)
      *Lxx *= h;
    if(Lp)
      *Lp *= h;
    if(Lpp)
      *Lpp *= h;
    if(Lpx)
      *Lpx *= h;
    L *= h;
  }
  */
  /*
  if(Lp)
    cout << "Lp: " << endl << *Lp << endl;
  if(Lpp)
    cout << "Lpp: " << endl << *Lpp << endl;
  if(Lpx)
    cout << "Lpx: " << endl << *Lpx << endl;
  if(Lx)
    cout << "Lx: " << endl << *Lx << endl;
  if(Lxx)
    cout << "Lxx: " << endl << *Lxx << endl;
  */

  return L;
}
