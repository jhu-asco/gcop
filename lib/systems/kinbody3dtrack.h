#ifndef GCOP_KINBODY3DTRACK_H
#define GCOP_KINBODY3DTRACK_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include "se3.h"
#include <iostream>
#include "utils.h"
#include "kinbody3d.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  
  /**
   * Pose-track in 2D. This is the simplese possible definition assuming point (unprojected) features
   * e.g. from stereo/laser scans with uniform spherical covariance. The framework can
   * be easily extended to any general 2d/3d non-uniform or projected features as well.
   */
  template <int _nu = 6>
  class Kinbody3dTrack {
  public:
    typedef Matrix<double, _nu, 1> Vectorud;
    typedef Matrix<double, 6, 1> Vector6d;
    /**
     * Pose graph (using simulated features around a circular track)
     * @param sys body3d system
     * @param nf number of features (landmarks)
     * @param t0 start time
     * @param tf end time
     * @param r radius
     * @param odometry is odometry available
     * @param extforce treat constant external force in x-y as a parameter
     * @param forces use uncertain forces in the cost function
     */
    Kinbody3dTrack(Kinbody3d<_nu> &sys, int nf, double vd0, double t0, double tf,
                double r = 25,
                bool extforce = false,
                bool forces = false);

    virtual void Get(Matrix4d &x, double vd, double t) const;

    virtual void MakeTrue();

    virtual void Add(const Vectorud &u, const Matrix4d &x, double h);

    /**
     * Add a new state/control to estimation vector
     * @param u commanded control 
     * @param x true state (to generate measurements)
     * @param h time-step 
     */
    virtual void Add2(const Vectorud &u, const Matrix4d &x, double h);

    Kinbody3d<_nu> &sys;     ///< system

    double t0;       ///< initial time around track
    double tf;       ///< final time around track
    double r;       ///< track radius
    double w;      ///< track width
    double h;      ///< track height
    double dmax;    ///< sensor radius

    bool extforce;         ///< is there a constant external force that need to be estimated

    bool forces;           ///< is there a constant external force that need to be estimated
    
    vector<double> ts;     ///< sequence of times (N+1 vector)

    vector<Matrix4d> xs;      ///< sequence of states x=(g,v) (N+1 vector)

    vector<Vectorud> us;   ///< sequence of inputs (N vector)
    
    vector<Vector3d> ls;   ///< landmark positions
    vector<bool> observed; ///< was it observed

    vector<Matrix4d> xos;     ///< unoptimized trajectory
    vector<Vectorud> uos;     ///< unoptimized trajectory

    bool init; 
    VectorXd p;            ///< observed landmark positions (2*m vector) and external force parameters (2 vector)    
    vector<int> pis;       ///< observed landmark indices (m vector)

    double pr;     ///< feature radius

    vector< vector<pair<int, Vector3d> > > Is;  ///< N+1-vector of vectors of visible pose-to-feature indices

    vector< vector<pair<int, Vector3d> > > Js;  ///< nf-vector of vectors of feature-to-pose indices indices
    vector<int> cis;       ///< if observed then this should map into the corresponding value in p

    double cp;             ///< noise covariance of poses (assume spherical)

    Vector6d cv;           ///< noise covariance of odometry/gyro (assume diagonal)

    Vectorud cw;           ///< noise covariance of control forces (assume diagonal)

  };

template <int _nu>
Kinbody3dTrack<_nu>::Kinbody3dTrack(Kinbody3d<_nu> &sys, int nf, double vd0, double t0, double tf,
                         double r,
                         bool extforce,
                         bool forces) : 
  sys(sys), t0(t0), tf(tf), r(r), w(4), h(4), dmax(20),
  extforce(extforce), forces(forces),
  ts(1,t0), ls(nf), observed(nf, false), p(extforce*3), pr(.75),
  Is(1), Js(nf), cis(nf), cp(.01)
{
  Matrix4d x;
  Get(x, vd0, t0);
  xs.push_back(x);
  xos.push_back(x);
  cw = MatrixXd::Constant(_nu, 1, 0.1);
}


template <int _nu>
void Kinbody3dTrack<_nu>::Add(const Vectorud &u, const Matrix4d &x, double h)
{
  // add noisy controls/state to list, but use true controls/state for generating measurements
  Eigen::Matrix4d xn;
  double t = ts.back();
  
  // noisy controls
  Vectorud un;
  for(int i = 0; i < un.size(); i++)
  {
    un[i] = u[i] + sqrt(cw[i])*random_normal();
  }

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
      Vector3d z = R.transpose()*(pf - px);
      z(0) += sqrt(cp)*random_normal();
      z(1) += sqrt(cp)*random_normal();
      z(2) += sqrt(cp)*random_normal();

      if (!observed[l]) {
        if (!init) {
          p.resize(3);
          init= true;
        } else {
          p.conservativeResize(p.size() + 3);
        }
        const Matrix3d &Rn = xs.back().block<3,3>(0,0);
        const Vector3d &pxn = xs.back().block<3,1>(0,3);
         
        // Initialize landmark position from estimated state
        p.tail<3>() = Rn*z + pxn;
        pis.push_back(l);
        observed[l] = true;
        cis[l] = p.size()/3-1;
      }

      //        zs[k].push_back();      // add feature l to pose k
      
      Is[k].push_back(make_pair(l, z)); // add feature l to pose k      
      Js[l].push_back(make_pair(k, z)); // add pose k to feature l
      
      //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
    }
  }
}


// add commanded control u, after which true noisy state x occured
template <int _nu>
void Kinbody3dTrack<_nu>::Add2(const Vectorud &u, const Matrix4d &x, double h)
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
      z(0) += sqrt(cp)*random_normal();
      z(1) += sqrt(cp)*random_normal();
      z(2) += sqrt(cp)*random_normal();
      
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
         
        // Initialize landmark position from estimated state
        p.tail<3>() = Rn*z + pxn;
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
template <int _nu>
void Kinbody3dTrack<_nu>::MakeTrue()
{
  //  srand (time(NULL));
  //  int nf = (p.size() - extforce*2)/2;
  int i = 0;  
  for (int l =0; l < ls.size();) {
    double a = RND*2*M_PI;
    if (a>M_PI && a<1.5*M_PI || a>0 && a<M_PI/2)
      continue;
    double z = (RND-0.5)*1.5*this->h;
    double r_rand = 3*w*(RND-0.5) + r; 
    ls[l] = Vector3d(r_rand*cos(a), r_rand*sin(a), z);
    ++l;
  }
}


// Get the desired state on the track for some time t
template <int _nu>
void Kinbody3dTrack<_nu>::Get(Matrix4d &x, double vd, double t) const
{
  double a = 1.3*t/tf*2*M_PI;
  x.setIdentity();
  x(0,0) = cos(a+M_PI/2.); x(0,1) = -sin(a+M_PI/2.);
  x(1,0) = sin(a+M_PI/2.); x(1,1) = cos(a+M_PI/2.);
 
  x(0,3) = r*cos(a); x(1,3) = r*sin(a); x(2,3) = 0;
}
}

#endif
