#include "kinbodyprojtrack.h"
#include "se3.h"
#include <iostream>
#include "utils.h"
#include <time.h>

using namespace gcop;
using namespace Eigen;
using namespace std;

KinbodyProjTrack::KinbodyProjTrack(Kinbody3d &sys, int nf, double vd0, double t0, double tf,
                         double r,
                         bool extforce,
                         bool forces) : 
  Kinbody3dTrack(sys, nf, vd0, t0, tf, r, extforce, forces)
{
  Matrix4d x;
  Get(x, vd0, t0);
  xs.clear();
  xos.clear();
  xs.push_back(x);
  xos.push_back(x);
}


void KinbodyProjTrack::Add(const Vector6d &u, const Matrix4d &x, double h)
{
  // add noisy controls/state to list, but use true controls/state for generating measurements
  Eigen::Matrix4d xn;
  double t = ts.back();
  
  // noisy controls
  Vector6d un;
  un << u[0] + sqrt(cw[0])*random_normal(),
              u[1] + sqrt(cw[1])*random_normal(),
              u[2] + sqrt(cw[2])*random_normal(),
              u[3] + sqrt(cw[3])*random_normal(),
              u[4] + sqrt(cw[4])*random_normal(),
              u[5] + sqrt(cw[5])*random_normal();


  sys.Step(xn, t, xs.back(), un, h);
  us.push_back(un);
  xs.push_back(xn);
  ts.push_back(t + h);

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
  
  const Matrix3d &R = x.block<3,3>(0,0);
  const Vector3d &px = x.block<3,1>(0,3);
  
  Is.resize(Is.size()+1);
  assert(k == Is.size()-1);

  for (int l = 0; l < ls.size(); ++l) {
    const Vector3d &pf = ls[l];
    //p.segment<2>(2*extforce + 2*l);
    
    double d = (px - pf).norm();
    if (d < dmax) {
      // add to feature vector if not observed
      if (!observed[l]) {
        if (!init) {
          p.resize(3);
          init= true;
        } else {
          p.conservativeResize(p.size() + 3);
        }
        p.tail<3>() = pf;
        pis.push_back(l);
        observed[l] = true;
        cis[l] = p.size()/3-1;
      }

      Vector3d z = R.transpose()*(pf - px);
      z /= z.norm();
      z(0) += sqrt(cp)*random_normal();
      z(1) += sqrt(cp)*random_normal();
      z(2) += sqrt(cp)*random_normal();
      z /= z.norm();
      //        zs[k].push_back();      // add feature l to pose k
      
      Is[k].push_back(make_pair(l, z)); // add feature l to pose k      
      Js[l].push_back(make_pair(k, z)); // add pose k to feature l
      
      //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
    }
  }
}


// add commanded control u, after which true noisy state x occured
void KinbodyProjTrack::Add2(const Vector6d &u, const Matrix4d &x, double h)
{
  // add noisy controls/state to list, but use true controls/state for generating measurements
  Eigen::Matrix4d xn;
  double t = ts.back();
  

  sys.Step(xn, t, xs.back(), u, h);
  us.push_back(u);
  xs.push_back(xn);
  ts.push_back(t + h);

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
  
  const Matrix3d &R = x.block<3,3>(0,0);
  const Vector3d &px = x.block<3,1>(0,3);
  
  Is.resize(Is.size()+1);
  assert(k == Is.size()-1);

  for (int l = 0; l < ls.size(); ++l) {
    const Vector3d &pf = ls[l];
    //p.segment<2>(2*extforce + 2*l);
    
    double d = (px - pf).norm();
    if (d < dmax) {

      Vector3d z = R.transpose()*(pf - px);
      z /= z.norm();
      z(0) += sqrt(cp)*random_normal();
      z(1) += sqrt(cp)*random_normal();
      z(2) += sqrt(cp)*random_normal();
      z /= z.norm();
      
      // add to feature vector if not observed
      if (!observed[l]) {
        if (!init) {
          p.resize(3);
          init= true;
        } else {
          p.conservativeResize(p.size() + 3);
        }
        const Matrix3d &Rn = xs.back().block<3,3>(0,0);
        const Vector3d &pxn = xs.back().block<3,1>(0,3);
         
        //cout << "feature z:" << endl << z << endl;
        //cout << "est x:" << endl << pxn << endl << Rn << endl;
        //cout << "true x:" << endl << px << endl << R;

        // Initialize landmark position from estimated state
        p.tail<3>() = Rn*dmax*z + pxn;
        // Initialize landmark position with its true position
        //p.tail<3>() = pf;
        pis.push_back(l);
        observed[l] = true;
        cis[l] = p.size()/3-1;
      }
      
      Is[k].push_back(make_pair(l, z)); // add feature l to pose k      
      Js[l].push_back(make_pair(k, z)); // add pose k to feature l
      
      //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
    }
  }
}


// Makes the ground truth features to be observed
void KinbodyProjTrack::MakeTrue()
{
  //  srand (time(NULL));
  //  int nf = (p.size() - extforce*2)/2;
  int i = 0;  
  int bound = 40;
  for (int l =0; l < ls.size();) {
    //double a = RND*bound;
    double a = RND*2*M_PI;
    if (a> 0.2 + M_PI && a<0.2 + 1.5*M_PI || a>0.2 && a<0.2 + 0.5*M_PI)
      continue;
    double z = (RND-0.5)*1.5*this->h;
    double r_rand = 3*w*(RND-0.5) + r; 
    ls[l] = Vector3d(r_rand*cos(a), r_rand*sin(a), z);
    //ls[l] = Vector3d(a-bound/2., 5+z, z);
    ++l;
  }
}


void KinbodyProjTrack::Optp(VectorXd &p, const vector<Matrix4d> &xs)
{
  int nf = (p.size() - extforce*3)/3;
  for (int l = 0; l < nf; ++l) {
    const vector< pair<int,Vector3d> > &J = Js[l];
    assert(J.size());
    Vector3d pf;// = Vector3d::Zero();
		pf.setZero();
    for (int j = 0; j < J.size(); ++j) {
      int k = J[j].first;
      const Vector3d &z = J[j].second;
      const Matrix3d &R = xs[k].block<3,3>(0,0);
      const Vector3d &x = xs[k].block<3,1>(0,3);
      pf = pf + x + R*z;
    }
    p.segment<3>(3*extforce + 3*l) = pf/J.size();
  }
}

// Get the desired state on the track for some time t
void KinbodyProjTrack::Get(Matrix4d &x, double vd, double t) const
{
  double a = 1.3*2*M_PI*t/tf;
  //double a = 20*t/tf;
  x.setIdentity();
  
  x(0,0) = cos(a+M_PI/2.); x(0,1) = -sin(a+M_PI/2.);
  x(1,0) = sin(a+M_PI/2.); x(1,1) = cos(a+M_PI/2.);
 
  x(0,3) = r*cos(a); x(1,3) = r*sin(a); x(2,3) = 0;
  
  //x(0,3) = 3*a;
}
