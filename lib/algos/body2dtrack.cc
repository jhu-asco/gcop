#include "body2dtrack.h"
#include "se2.h"
#include <iostream>
#include "utils.h"
#include <time.h>

using namespace gcop;
using namespace Eigen;
using namespace std;

Body2dTrack::Body2dTrack(Body2d &sys, int nf, double t0, double tf,
                         double r,
                         bool odometry,
                         bool extforce,
                         bool forces) : 
  sys(sys), t0(t0), tf(tf), r(r), w(4), dmax(10),
  odometry(odometry), extforce(extforce), forces(forces),
  ts(1,t0), ls(nf), observed(nf, false), p(extforce*2), pr(.75),
  Is(1), Js(nf), cis(nf), cp(.01), cv(0.05, 0.1, 0.1), cw(.5, 2, 2)
{
  M3V3d x;
  Get(x, 5, t0);
  xs.push_back(x);
  xos.push_back(x);
  vs.push_back(x.second);
}


void Body2dTrack::Add(const Vector3d &u, const M3V3d &x, double h)
{
  // add noisy controls/state to list, but use true controls/state for generating measurements
  M3V3d xn;
  double t = ts.back();
  
  // noisy controls
  Vector3d un(u[0] + sqrt(cw[0])*random_normal(),
              u[1] + sqrt(cw[1])*random_normal(),
              u[2] + sqrt(cw[2])*random_normal());

  sys.Step(xn, t, xs.back(), un, h);
  us.push_back(un);
  xs.push_back(xn);
  ts.push_back(t + h);

  if (odometry) {
    Vector3d vn(sqrt(cv[0])*random_normal(),
                sqrt(cv[1])*random_normal(),
                sqrt(cv[2])*random_normal());
    
    vs.push_back(x.second + vn); // set noisy gyro measurements
  }

  sys.Step(xn, t, xos.back(), un, h);
  xos.push_back(xn);
  uos.push_back(u);
  
  // time-stage 
  int k = xs.size()-1;

  /*
  sys.force->fext << 0, .5, 0;
  if (extforce) {
    p.head<2>() << .5, 0;
  }
  */
  
  /* // add noise
  for (int k = 0; k < N; ++k) {
    double t = k*h;
    Get(xs[k+1], r, vd, t+h);
    Vector3d f;
    us[k] = sys.I.cwiseProduct(xs[k+1].second - xs[k].second)/h + sys.force->D.cwiseProduct(xs[k].second);
    // add the constant external force
    us[k].tail<2>() -= xs[k].first.topLeftCorner<2,2>().transpose()*sys.force->fext.tail<2>();
  }
  */

  //  srand (time(NULL));
  
  const Matrix2d &R = x.first.topLeftCorner<2,2>();
  const Vector2d &px = x.first.block<2,1>(0,2);
  
  Is.resize(Is.size()+1);
  assert(k == Is.size()-1);

  for (int l = 0; l < ls.size(); ++l) {
    const Vector2d &pf = ls[l];
    //p.segment<2>(2*extforce + 2*l);
    
    double d = (px - pf).norm();
    if (d < dmax) {
      // add to feature vector if not observed
      if (!observed[l]) {
        if (!init) {
          p.resize(2);
          init= true;
        } else {
          p.conservativeResize(p.size() + 2);
        }
        p.tail<2>() = pf;
        pis.push_back(l);
        observed[l] = true;
        cis[l] = p.size()/2-1;
      }

      Vector2d z = R.transpose()*(pf - px);
      z(0) += sqrt(cp)*random_normal();
      z(1) += sqrt(cp)*random_normal();
      
      //        zs[k].push_back();      // add feature l to pose k
      
      Is[k].push_back(make_pair(l, z)); // add feature l to pose k      
      Js[l].push_back(make_pair(k, z)); // add pose k to feature l
      
      //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
    }
  }
}


// add commanded control u, after which true noisy state x occured
void Body2dTrack::Add2(const Vector3d &u, const M3V3d &x, double h)
{
  // add noisy controls/state to list, but use true controls/state for generating measurements
  M3V3d xn;
  double t = ts.back();
  
  sys.Step(xn, t, xs.back(), u, h);
  us.push_back(u);
  xs.push_back(xn);
  ts.push_back(t + h);

  if (odometry) {
    Vector3d vn(sqrt(cv[0])*random_normal(),
                sqrt(cv[1])*random_normal(),
                sqrt(cv[2])*random_normal());
    
    vs.push_back(x.second + vn); // set noisy gyro measurements
  }

  sys.Step(xn, t, xos.back(), u, h);
  xos.push_back(xn);
  uos.push_back(u);
  
  // time-stage 
  int k = xs.size()-1;

  /*
  sys.force->fext << 0, .5, 0;
  if (extforce) {
    p.head<2>() << .5, 0;
  }
  */
  
  /* // add noise
  for (int k = 0; k < N; ++k) {
    double t = k*h;
    Get(xs[k+1], r, vd, t+h);
    Vector3d f;
    us[k] = sys.I.cwiseProduct(xs[k+1].second - xs[k].second)/h + sys.force->D.cwiseProduct(xs[k].second);
    // add the constant external force
    us[k].tail<2>() -= xs[k].first.topLeftCorner<2,2>().transpose()*sys.force->fext.tail<2>();
  }
  */

  //  srand (time(NULL));
  
  const Matrix2d &R = x.first.topLeftCorner<2,2>();
  const Vector2d &px = x.first.block<2,1>(0,2);
  
  Is.resize(Is.size()+1);
  assert(k == Is.size()-1);

  for (int l = 0; l < ls.size(); ++l) {
    const Vector2d &pf = ls[l];
    //p.segment<2>(2*extforce + 2*l);
    
    double d = (px - pf).norm();
    if (d < dmax) {
      // add to feature vector if not observed
      if (!observed[l]) {
        if (!init) {
          p.resize(2);
          init= true;
        } else {
          p.conservativeResize(p.size() + 2);
        }
        p.tail<2>() = pf;
        pis.push_back(l);
        observed[l] = true;
        cis[l] = p.size()/2-1;
      }

      Vector2d z = R.transpose()*(pf - px);
      z(0) += sqrt(cp)*random_normal();
      z(1) += sqrt(cp)*random_normal();
      
      //        zs[k].push_back();      // add feature l to pose k
      
      Is[k].push_back(make_pair(l, z)); // add feature l to pose k      
      Js[l].push_back(make_pair(k, z)); // add pose k to feature l
      
      //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
    }
  }
}


void Body2dTrack::MakeTrue()
{
  //  srand (time(NULL));
  //  int nf = (p.size() - extforce*2)/2;
  int i = 0;  
  for (int l =0; l < ls.size();) {
    double a = RND*2*M_PI;
    if (a>M_PI && a<1.5*M_PI || a>0 && a<M_PI/2)
      continue;
    ls[l] = 1.5*r*RND*Vector2d(cos(a), sin(a));
    ++l;
  }
}


void Body2dTrack::Optp(VectorXd &p, const vector<M3V3d> &xs)
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


void Body2dTrack::Get(M3V3d &x, double vd, double t) const
{
  double a = 1.3*t/tf*2*M_PI;
  x.first.setIdentity();
  SE2::Instance().q2g(x.first, Vector3d(a + M_PI/2, r*cos(a), r*sin(a)));
  x.second << vd/r, vd, 0;
}
