#include "posegraph2d.h"
#include "se2.h"
#include <iostream>
#include "utils.h"
#include <time.h>

using namespace gcop;
using namespace Eigen;
using namespace std;

Posegraph2d::Posegraph2d(int N, int nf) : 
  ts(N+1), gs(N+1), us(N), p(2*nf), Is(N+1), Js(nf), cp(.01)
{
}

void Posegraph2d::Synthesize(Posegraph2d &pgt, Posegraph2d &pgn, double tf)
{
  int N = pgt.us.size();
  int nf = pgt.p.size()/2;
  double h = tf/N;

  pgt.gs[0].setIdentity();
  pgn.gs[0].setIdentity();
  
  pgt.ts[0] = 0;
  pgn.ts[0] = 0;

  Matrix3d m;
  
  for (int k = 0; k < N; ++k) {
    double t = k*h;

    pgt.us[k] << cos(2*t), 1, 0;
    //    pgt.us[k] << 0, 1, 0;
    SE2::Instance().cay(m, h*pgt.us[k]); 
    pgt.gs[k+1] = pgt.gs[k]*m;
    pgt.ts[k+1] = (k+1)*h;
    
    //    pgn.us[k] = pgt.us[k] + Vector3d(0.3*random_normal(), random_normal(), 0); // add drift
    pgn.us[k] = pgt.us[k] + Vector3d(-.1, 0, 0); // add drift

    SE2::Instance().cay(m, h*pgn.us[k]); 
    pgn.gs[k+1] = pgn.gs[k]*m;    
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
    const Matrix2d &R = pgt.gs[k].topLeftCorner<2,2>();
    const Vector2d &x = pgt.gs[k].block<2,1>(0,2);
    
    for (int j = 0; j < nfp; ++j) {
      double a = 2*M_PI*RND;
      pgt.p.segment<2>(2*l) = x + dmax*RND*Vector2d(cos(a), sin(a));
      ++l;
    }
  }
  
  srand (time(NULL));
 
  for (int k = 0; k <= N; ++k) {
    const Matrix2d &R = pgt.gs[k].topLeftCorner<2,2>();
    const Vector2d &x = pgt.gs[k].block<2,1>(0,2);

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

  pgn.Optp(pgn.p, pgn.gs);
}


void Posegraph2d::Synthesize2(Posegraph2d &pgt, Posegraph2d &pgn, double tf)
{
  int N = pgt.us.size();
  int nf = pgt.p.size()/2;
  double h = tf/N;

  pgt.gs[0].setIdentity();
  pgn.gs[0].setIdentity();
  
  pgt.ts[0] = 0;
  pgn.ts[0] = 0;

  Matrix3d m;
  
  for (int k = 0; k < N; ++k) {
    double t = k*h;

    //    pgt.us[k] << cos(2*t), 1, 0;
    double w= .5;
    if (t*w < 2*M_PI)
      pgt.us[k] << w, 1, 0;
    else
      pgt.us[k] << -w, 1, 0;

    
    //    pgt.us[k] << 0, 1, 0;
    SE2::Instance().cay(m, h*pgt.us[k]); 
    pgt.gs[k+1] = pgt.gs[k]*m;
    pgt.ts[k+1] = (k+1)*h;
    
    pgn.us[k] = pgt.us[k] + Vector3d(-0.1, 2*random_normal(), 2*random_normal()); // add drift
    SE2::Instance().cay(m, h*pgn.us[k]); 
    pgn.gs[k+1] = pgn.gs[k]*m;    
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
    const Matrix2d &R = pgt.gs[k].topLeftCorner<2,2>();
    const Vector2d &x = pgt.gs[k].block<2,1>(0,2);
    
    for (int j = 0; j < nfp; ++j) {
      double a = 2*M_PI*RND;
      pgt.p.segment<2>(2*l) = x + dmax*RND*Vector2d(cos(a), sin(a));
      ++l;
    }
  }
  
  srand (time(NULL));
 
  for (int k = 0; k <= N; ++k) {
    const Matrix2d &R = pgt.gs[k].topLeftCorner<2,2>();
    const Vector2d &x = pgt.gs[k].block<2,1>(0,2);

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

  pgn.Optp(pgn.p, pgn.gs);
}


void Posegraph2d::Optp(VectorXd &p, const vector<Matrix3d> &gs)
{
  int nf = p.size()/2;
  for (int l = 0; l < nf; ++l) {
    const vector< pair<int,Vector2d> > &J = Js[l];
    assert(J.size());
    Vector2d pf = Vector2d::Zero();
    for (int j = 0; j < J.size(); ++j) {
      int k = J[j].first;
      const Vector2d &z = J[j].second;
      const Matrix2d &R = gs[k].topLeftCorner<2,2>();
      const Vector2d &x = gs[k].block<2,1>(0,2);
      pf = pf + x + R*z;
    }
    p.segment<2>(2*l) = pf/J.size();
  }
}
