#include "body2dgraph.h"
#include "se2.h"
#include <iostream>
#include "utils.h"
#include <time.h>

using namespace gcop;
using namespace Eigen;
using namespace std;

Body2dGraph::Body2dGraph(Body2d<> &sys, int N, int nf, 
                         bool odometry, 
                         bool extforce,
                         bool forces) : 
  sys(sys), odometry(odometry), extforce(extforce), forces(forces), fext(0,0,0),
  ts(N+1), xs(N+1,make_pair(Matrix3d::Identity(), Vector3d::Zero())), us(N), 
  p(extforce*2 + 2*nf), 
  Is(N+1), Js(nf), vs(N+1), cp(.01), cv(0.01, 0.1, 0.1), cw(.005, .03, .03)
{
}

void Body2dGraph::Synthesize(Body2dGraph &pgt, Body2dGraph &pgn, double tf)
{
  int N = pgt.us.size();
  int nf = pgt.p.size()/2;
  double h = tf/N;

  pgt.xs[0].first.setIdentity();
  pgt.xs[0].second.setZero();
  pgn.xs[0].first.setIdentity();
  pgn.xs[0].second.setZero();
  
  pgt.ts[0] = 0;
  pgn.ts[0] = 0;

  Matrix3d m;
  
  for (int k = 0; k < N; ++k) {
    double t = k*h;

    pgt.xs[k+1].second << cos(2*t), 1, 0;
    //    pgt.us[k] << 0, 1, 0;
    SE2::Instance().cay(m, h*pgt.xs[k+1].second); 
    pgt.xs[k+1].first = pgt.xs[k].first*m;
    pgt.us[k] = (pgt.xs[k+1].second - pgt.xs[k].second)/h;
    pgt.ts[k+1] = (k+1)*h;
    
    //    pgn.us[k] = pgt.us[k] + Vector3d(0.3*random_normal(), random_normal(), 0); // add drift
    pgn.xs[k+1].second = pgt.xs[k+1].second + Vector3d(-.1, 0, 0); // add drift

    SE2::Instance().cay(m, h*pgn.xs[k+1].second); 
    pgn.xs[k+1].first = pgn.xs[k].first*m;    
    pgn.us[k] = (pgn.xs[k+1].second - pgn.xs[k].second)/h;
    pgn.ts[k+1] = (k+1)*h;
  }
  
  int nfp = nf/(N+1); // features per pose

  int l = 0;        // current feature index

  double dmax = .75;  // visibility radius

  // create the features
  for (int k = 0; k <= N; ++k) {
    if (l >= nf)
      break;
    
    //    cout << "k=" << k << " " << gs[k] << endl;
    const Matrix2d &R = pgt.xs[k].second.topLeftCorner<2,2>();
    const Vector2d &x = pgt.xs[k].second.block<2,1>(0,2);
    
    for (int j = 0; j < nfp; ++j) {
      double a = 2*M_PI*RND;
      pgt.p.segment<2>(2*l) = x + dmax*RND*Vector2d(cos(a), sin(a));
      ++l;
    }
  }
  
  srand (time(NULL));
 
  for (int k = 0; k <= N; ++k) {
    const Matrix2d &R = pgt.xs[k].second.topLeftCorner<2,2>();
    const Vector2d &x = pgt.xs[k].second.block<2,1>(0,2);

    for (int l = 0; l < nf; ++l) {
      const Vector2d &pf = pgt.p.segment<2>(2*l);

      double d = (x - pf).norm();
      if (d < dmax) {
        Vector2d z = R.transpose()*(pf - x);
        z(0) += sqrt(pgt.cp)*random_normal();
        z(1) += sqrt(pgt.cp)*random_normal();

        //        zs[k].push_back();      // add feature l to pose k
        
        pgt.Is[k].push_back(make_pair(l, z)); // add feature l to pose k

        pgt.Js[l].push_back(make_pair(k, z)); // add pose k to feature l

        //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
      }
    }
  }

  /*
  cout << "N=" << N << " " << gs.size() << endl;
  cout << p.size() << endl;
  cout << "p=" << p << endl;
  cout << "gs[N]=" << gs[N];
  cout << "Is[N]=" << Is[N][0].first << " " << Is[N][0].second << endl;
 
  for (int k=0; k< Js[Is[N][0].first].size(); ++k)
    cout << "Js[x]=" << Js[Is[N][0].first][k].first << " " << Js[Is[N][0].first][k].second << endl;
  */
  //  exit(0);
  
  // estimate noisy features

  pgn.Is = pgt.Is;
  pgn.Js = pgt.Js;

  pgn.Optp(pgn.p, pgn.xs);
}


void Body2dGraph::Synthesize2(Body2dGraph &pgt, Body2dGraph &pgn, double tf)
{
  //  srand (time(NULL));

  int N = pgt.us.size();
  int nf = pgt.p.size()/2;
  double h = tf/N;

  pgt.xs[0].first.setIdentity();
  pgt.xs[0].second.setZero();
  pgn.xs[0].first.setIdentity();
  pgn.xs[0].second.setZero();
  
  pgt.ts[0] = 0;
  pgn.ts[0] = 0;

  pgn.vs[0] = pgn.xs[0].second; // set gyro measurements

  Matrix3d m;

  if (pgt.extforce) {
    pgt.p.tail(2) << -.1, 0;
    pgn.p.tail(2) << 0, 0;
  }
  
  for (int k = 0; k < N; ++k) {
    double t = k*h;

    //    pgt.us[k] << cos(2*t), 1, 0;
    double w= .5;
    if (t*w < 2*M_PI)
      pgt.xs[k+1].second << w, 1, 0;
    else
      pgt.xs[k+1].second << -w, 1, 0;
    
    //    pgt.us[k] << 0, 1, 0;
    SE2::Instance().cay(m, h*pgt.xs[k+1].second); 
    pgt.xs[k+1].first = pgt.xs[k].first*m;
    pgt.us[k] = (pgt.xs[k+1].second - pgt.xs[k].second)/h + pgt.sys.force->D.cwiseProduct(pgt.xs[k].second);
    pgt.ts[k+1] = (k+1)*h;
    
    pgn.xs[k+1].second = pgt.xs[k+1].second + Vector3d(-0.05, random_normal(), random_normal()); // add drift

    //    pgn.xs[k+1].second = pgt.xs[k+1].second + Vector3d(sqrt(pgn.cv[0])*random_normal(), sqrt(pgn.cv[1])*random_normal(), sqrt(pgn.cv[2])*random_normal());

    SE2::Instance().cay(m, h*pgn.xs[k+1].second); 
    pgn.xs[k+1].first = pgn.xs[k].first*m;    
    pgn.us[k] = (pgn.xs[k+1].second - pgn.xs[k].second)/h + pgn.sys.force->D.cwiseProduct(pgn.xs[k].second);
    pgn.ts[k+1] = (k+1)*h;

    pgn.vs[k+1] = pgn.xs[k+1].second; // set gyro measurements
  }
  
  int nfp = nf/(N+1); // features per pose

  int l = 0;        // current feature index

  double dmax = .75;  // visibility radius

  // create the features
  for (int k = 0; k <= N; ++k) {
    if (l >= nf)
      break;
    
    //    cout << "k=" << k << " " << gs[k] << endl;
    const Matrix2d &R = pgt.xs[k].first.topLeftCorner<2,2>();
    const Vector2d &x = pgt.xs[k].first.block<2,1>(0,2);
    
    for (int j = 0; j < nfp; ++j) {
      double a = 2*M_PI*RND;
      pgt.p.segment<2>(2*l) = x + dmax*RND*Vector2d(cos(a), sin(a));
      ++l;
    }
  }
  
  srand (time(NULL));
 
  for (int k = 0; k <= N; ++k) {
    const Matrix2d &R = pgt.xs[k].first.topLeftCorner<2,2>();
    const Vector2d &x = pgt.xs[k].first.block<2,1>(0,2);

    for (int l = 0; l < nf; ++l) {
      const Vector2d &pf = pgt.p.segment<2>(2*l);

      double d = (x - pf).norm();
      if (d < dmax) {
        Vector2d z = R.transpose()*(pf - x);
        z(0) += sqrt(pgt.cp)*random_normal();
        z(1) += sqrt(pgt.cp)*random_normal();

        //        zs[k].push_back();      // add feature l to pose k
        
        pgt.Is[k].push_back(make_pair(l, z)); // add feature l to pose k

        pgt.Js[l].push_back(make_pair(k, z)); // add pose k to feature l

        //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
      }
    }
  }

  pgn.Is = pgt.Is;
  pgn.Js = pgt.Js;

  pgn.Optp(pgn.p, pgn.xs);
}


void Body2dGraph::Synthesize3(Body2dGraph &pgt, Body2dGraph &pgn, double tf)
{
  srand (time(NULL));

  int N = pgt.us.size();
  int nf = (pgt.p.size() - pgt.extforce*2)/2;
  double h = tf/N;

  pgt.ts[0] = 0;
  pgt.xs[0].first.setIdentity();
  pgt.xs[0].second.setZero();

  pgn.ts[0] = 0;
  pgn.xs[0].first.setIdentity();
  pgn.xs[0].second.setZero();  

  pgn.vs[0] = pgn.xs[0].second; // set gyro measurements

  Matrix3d m;

  // pgt.sys.force->fext << 0, -1, 0;
  if (pgt.extforce) {
    pgt.p.head<2>() << .5, 0;
    pgn.p.head<2>() << .25, 0;
  }
  
  for (int k = 0; k < N; ++k) {
    double t = k*h;

    // TRUE TRAJECTORY
    double w = .5;
    if (t*w < 2*M_PI)
      pgt.xs[k+1].second << w, 1, 0;
    else
      pgt.xs[k+1].second << -w, 1, 0;
    
    //    pgt.us[k] << 0, 1, 0;
    SE2::Instance().cay(m, h*pgt.xs[k+1].second);
    pgt.ts[k+1] = (k+1)*h;
    pgt.xs[k+1].first = pgt.xs[k].first*m;
    Vector3d f;
    pgt.us[k] = pgt.sys.I.cwiseProduct(pgt.xs[k+1].second - pgt.xs[k].second)/h + pgt.sys.force->D.cwiseProduct(pgt.xs[k].second);
    // add the unknown constant external force
    pgt.us[k][0] -= pgt.fext[0]; 
    pgt.us[k].tail<2>() -= pgt.xs[k].first.topLeftCorner<2,2>().transpose()*pgt.fext.tail<2>();
    
    // NOISY TRAJECTORY
    pgn.ts[k+1] = (k+1)*h;

    // control forces noise
    Vector3d un(sqrt(pgt.cw[0])*random_normal(),
                sqrt(pgt.cw[1])*random_normal(),
                sqrt(pgt.cw[2])*random_normal());
    
    pgn.us[k] = pgt.us[k] + un;
    pgn.sys.Step(pgn.xs[k+1], t, pgn.xs[k], pgn.us[k], h);

    Vector3d vn(sqrt(pgt.cv[0])*random_normal(),
                sqrt(pgt.cv[1])*random_normal(),
                sqrt(pgt.cv[2])*random_normal());

    pgn.vs[k+1] = pgt.xs[k+1].second + vn; // set noisy gyro measurements
  }
  
  int nfp = nf/(N+1); // features per pose

  int l = 0;        // current feature index

  double dmax = .75;  // visibility radius

  // create the features
  for (int k = 0; k <= N; ++k) {
    if (l >= nf)
      break;
    
    //    cout << "k=" << k << " " << gs[k] << endl;
    const Matrix2d &R = pgt.xs[k].first.topLeftCorner<2,2>();
    const Vector2d &x = pgt.xs[k].first.block<2,1>(0,2);
    
    for (int j = 0; j < nfp; ++j) {
      double a = 2*M_PI*RND;
      pgt.p.segment<2>(2*pgt.extforce + 2*l) = x + dmax*RND*Vector2d(cos(a), sin(a));
      ++l;
    }
  }
  
  srand (time(NULL));
 
  for (int k = 0; k <= N; ++k) {
    const Matrix2d &R = pgt.xs[k].first.topLeftCorner<2,2>();
    const Vector2d &x = pgt.xs[k].first.block<2,1>(0,2);

    for (int l = 0; l < nf; ++l) {
      const Vector2d &pf = pgt.p.segment<2>(2*pgt.extforce + 2*l);

      double d = (x - pf).norm();
      if (d < dmax) {
        Vector2d z = R.transpose()*(pf - x);
        z(0) += sqrt(pgt.cp)*random_normal();
        z(1) += sqrt(pgt.cp)*random_normal();

        //        zs[k].push_back();      // add feature l to pose k
        
        pgt.Is[k].push_back(make_pair(l, z)); // add feature l to pose k

        pgt.Js[l].push_back(make_pair(k, z)); // add pose k to feature l

        //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
      }
    }
  }

  pgn.Is = pgt.Is;
  pgn.Js = pgt.Js;

  pgn.Optp(pgn.p, pgn.xs);
}


void Body2dGraph::Optp(VectorXd &p, const vector<Body2dState> &xs)
{
  int nf = (p.size() - extforce*2)/2;
  for (int l = 0; l < nf; ++l) {
    const vector< pair<int,Vector2d> > &J = Js[l];
    assert(J.size());
    Vector2d pf;// = Vector2d::Zero();
		pf.setZero();
    for (int j = 0; j < J.size(); ++j) {
      int k = J[j].first;
      const Vector2d &z = J[j].second;
      const Matrix2d &R = xs[k].first.topLeftCorner<2,2>();
      const Vector2d &x = xs[k].first.block<2,1>(0,2);
      pf = pf + x + R*z;
    }
    p.segment<2>(2*extforce + 2*l) = pf/J.size();
  }
}
