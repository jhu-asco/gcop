#include <assert.h>
#include "body2dtrackcost.h"
#include "body2dmanifold.h"
#include "rn.h"

using namespace std;
using namespace gcop;
using namespace Eigen;

Body2dTrackCost::Body2dTrackCost(double tf, const Body2dTrack &pg) :
  Cost(Body2dManifold::Instance(), Rn<3>::Instance(), tf), pg(pg)
{
  
}


static RowVector2d r2hat(const Vector2d &a)
{
  return RowVector2d(-a[1], a[0]);
}


static double cross2(Vector2d a, Vector2d b)
{
  return a[0]*b[1] - a[1]*b[0];
}


double Body2dTrackCost::Lp(double t, const M3V3d &x, const Vector3d &u, const VectorXd &p,
                          Vector6d *Lx, Matrix6d *Lxx,
                          Vector3d *Lu, Matrix3d *Luu,
                          Matrix63d *Lxu,
                          VectorXd *Lp, MatrixXd *Lpp, MatrixX6d *Lpx)
{
  double L = 0;

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
  
  const Matrix3d &g = x.first;

  const Matrix2d &R = g.topLeftCorner<2,2>(); // orientation
  const Vector2d &xp = g.block<2,1>(0,2);     // position

  const Vector3d &v = x.second;               // velocity
  
  int N = pg.Is.size() - 1;
  double h = this->tf/N;

  int k = (int)round(t/h);
  assert(k >=0 && k <= N);
  const vector< pair<int, Vector2d> > &I = pg.Is[k];

  //  cout << "Body2dtrackCost: k=" << k << " " << I.size() << endl;

  int i0 = 2*pg.extforce;
  for (int i = 0; i < I.size(); ++i) {
    int l = pg.cis[I[i].first];  // feature index
    const Vector2d &z = I[i].second;

    const Vector2d &pf = p.segment<2>(i0 + 2*l);
    Vector2d y = R.transpose()*(pf - xp);
    Vector2d r = (y - z)/pg.cp;

    L += r.dot(y - z)/2;         // feature cost

    if (Lx) {
      (*Lx)[0] = (*Lx)[0] - cross2(y, r);
      Lx->segment<2>(1) = Lx->segment<2>(1) - r;
    }

    if (Lxx) {
      (*Lxx)(0,0) += y.dot(y)/pg.cp;
      Lxx->block<1,2>(0,1) += (-r2hat(r) + r2hat(y)/pg.cp);
      //      Lxx->topRightCorner<1,2>() += (-r2hat(r) + r2hat(y)/pg.cp);
      Lxx->block<2,1>(1,0) += (-r2hat(r) + r2hat(y)/pg.cp).transpose();      
      //      Lxx->bottomLeftCorner<2,1>() += (-r2hat(r) + r2hat(y)/pg.cp).transpose();      
      Lxx->block<2,2>(1,1) += Matrix2d::Identity()/pg.cp;
      //      Lxx->bottomRightCorner<2,2>() += Matrix2d::Identity()/pg.cp;
    }

    if (Lp)
      Lp->segment<2>(i0 + 2*l) = R*r;

    if (Lpp)
      Lpp->block<2,2>(i0 + 2*l, i0 + 2*l) = Matrix2d::Identity()/pg.cp;

    if (Lpx) {
      Lpx->block<2,1>(i0 + 2*l, 0) = -R*r2hat(z).transpose()/pg.cp;
      Lpx->block<2,2>(i0 + 2*l, 1) = -R/pg.cp;
    }
  }

  if (pg.odometry) {           // odometry/gyro cost
    Vector3d Cdv;
    Vector3d dv = v - pg.vs[k];
    Cdv = dv.cwiseQuotient(pg.cv);
    L += Cdv.dot(dv)/2;
    
    if (Lx)
      Lx->segment<3>(3) += Cdv;

    if (Lxx)
      Lxx->block<3,3>(3,3) += Vector3d::Ones().cwiseQuotient(pg.cv).asDiagonal();    
  }

  if (pg.forces && k < N) {
    Vector3d du = u - pg.uos[k];
    Vector3d wdu = du.cwiseQuotient(pg.cw);

    L += du.dot(wdu)/2;
    if (Lu)
      *Lu = wdu;
    if (Luu)
      *Luu = Vector3d::Ones().cwiseQuotient(pg.cw).asDiagonal();
  }

  return L;
}
