#include <assert.h>
#include "bacost.h"
#include "kinbody2dmanifold.h"
#include "rn.h"

using namespace std;
using namespace gcop;
using namespace Eigen;

BaCost::BaCost(Kinbody2d &sys, double tf, const Posegraph2d &pg) : 
  Cost(sys, tf), pg(pg)
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


double BaCost::L(double t, const Matrix3d &g, const Vector3d &u, double h,
                 const VectorXd *p,
                  Vector3d *Lx, Matrix3d *Lxx,
                  Vector3d *Lu, Matrix3d *Luu,
                  Matrix3d *Lxu,
                  VectorXd *Lp, MatrixXd *Lpp, MatrixX3d *Lpx)
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
  
  const Matrix2d &R = g.topLeftCorner<2,2>();
  const Vector2d &x = g.block<2,1>(0,2);
  
  int N = pg.Is.size() - 1;
  //  double h = this->tf/N;

  int k = (h < 1e-16 ? N : (int)round(t/h));
  //  cout <<"h="<<h<<endl;
  assert(k >=0 && k <= N);
  const vector< pair<int, Vector2d> > &I = pg.Is[k];

  //  cout << "BaCost: k=" << k << " " << I.size() << endl;

  for (int i = 0; i < I.size(); ++i) {
    int l = I[i].first;
    const Vector2d &z = I[i].second;

    const Vector2d &pf = p->segment<2>(2*l);
    Vector2d y = R.transpose()*(pf - x);
    Vector2d r = (y - z)/pg.cp;

    L += r.dot(y - z)/2;         // cost

    if (Lx) {
      (*Lx)[0] = (*Lx)[0] - cross2(y, r);
      Lx->tail(2) = Lx->tail(2) - r;
    }

    if (Lxx) {
      (*Lxx)(0,0) += y.dot(y)/pg.cp;
      Lxx->topRightCorner<1,2>() += (-r2hat(r) + r2hat(y)/pg.cp);
      Lxx->bottomLeftCorner<2,1>() += (-r2hat(r) + r2hat(y)/pg.cp).transpose();      
      Lxx->bottomRightCorner<2,2>() += Matrix2d::Identity()/pg.cp;
    }

    if (Lp)
      Lp->segment<2>(2*l) = R*r;

    if (Lpp)
      Lpp->block<2,2>(2*l, 2*l) = Matrix2d::Identity()/pg.cp;

    if (Lpx) {
      Lpx->block<2,1>(2*l, 0) = -R*r2hat(z).transpose()/pg.cp;
      Lpx->block<2,2>(2*l, 1) = -R/pg.cp;    
    }
  }
  return L;
}

