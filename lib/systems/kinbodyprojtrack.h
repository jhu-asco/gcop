#ifndef GCOP_KINBODYPROJTRACK_H
#define GCOP_KINBODYPROJTRACK_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include "kinbody3d.h"
#include "kinbody3dtrack.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  
  /**
   * Pose-track in 2D. This is the simplese possible definition assuming point (unprojected) features
   * e.g. from stereo/laser scans with uniform spherical covariance. The framework can
   * be easily extended to any general 2d/3d non-uniform or projected features as well.
   */
  template <int _nu = 6>
  class KinbodyProjTrack : public Kinbody3dTrack<_nu> {
  public:
    typedef Matrix<double, _nu, 1> Vectorud;
    typedef Matrix<double, 6, _nu> Matrix6ud;
    typedef Matrix<double, 6, 6> Matrix6d;
    
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
    KinbodyProjTrack(Kinbody3d<_nu> &sys, int nf, double vd0, double t0, double tf,
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
  };

template <int _nu>
KinbodyProjTrack<_nu>::KinbodyProjTrack(Kinbody3d<_nu> &sys, int nf, double vd0, double t0, double tf,
                         double r,
                         bool extforce,
                         bool forces) : 
  Kinbody3dTrack<_nu>(sys, nf, vd0, t0, tf, r, extforce, forces)
{
  Matrix4d x;
  Get(x, vd0, t0);
  this->xs.clear();
  this->xos.clear();
  this->xs.push_back(x);
  this->xos.push_back(x);
}


template <int _nu>
void KinbodyProjTrack<_nu>::Add(const Vectorud &u, const Matrix4d &x, double h)
{
  // add noisy controls/state to list, but use true controls/state for generating measurements
  Eigen::Matrix4d xn;
  double t = this->ts.back();
  
  Vectorud un;
  // noisy controls
  for(int i = 0; i < u.size(); i++)
  {
    un[i] = u[i] + sqrt(this->cw[i])*random_normal();
  }

  this->sys.Step(xn, t, this->xs.back(), un, h);
  this->us.push_back(un);
  this->xs.push_back(xn);
  this->ts.push_back(t + h);

  this->sys.Step(xn, t, this->xos.back(), un, h);
  this->xos.push_back(xn);
  this->uos.push_back(u);
  
  // time-stage 
  int k = this->xs.size()-1;

  const Matrix3d &R = x.block<3,3>(0,0);
  const Vector3d &px = x.block<3,1>(0,3);
  
  this->Is.resize(this->Is.size()+1);
  assert(k == this->Is.size()-1);

  for (int l = 0; l < this->ls.size(); ++l) {
    const Vector3d &pf = this->ls[l];
    //p.segment<2>(2*extforce + 2*l);
    
    double d = (px - pf).norm();
    if (d < this->dmax) {
      Vector3d z = R.transpose()*(pf - px);
      z /= z.norm();
      z(0) += sqrt(this->cp)*random_normal();
      z(1) += sqrt(this->cp)*random_normal();
      z(2) += sqrt(this->cp)*random_normal();
      z /= z.norm();

      // add to feature vector if not observed
      if (!this->observed[l]) {
        if (!this->init) {
          this->p.resize(3);
          this->init= true;
        } else {
          this->p.conservativeResize(this->p.size() + 3);
        }
        const Matrix3d &Rn = this->xs.back().block(0,0,3,3);
        const Vector3d &pxn = this->xs.back().block(0,3,3,1);
         
        // Initialize landmark position from estimated state
        this->p.tail(3) = Rn*this->dmax*z + pxn;
        this->pis.push_back(l);
        this->observed[l] = true;
        this->cis[l] = this->p.size()/3-1;
      }

      //        zs[k].push_back();      // add feature l to pose k
      
      this->Is[k].push_back(make_pair(l, z)); // add feature l to pose k      
      this->Js[l].push_back(make_pair(k, z)); // add pose k to feature l
      
      //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
    }
  }
}


// add commanded control u, after which true noisy state x occured
template <int _nu>
void KinbodyProjTrack<_nu>::Add2(const Vectorud &u, const Matrix4d &x, double h)
{
  // add noisy controls/state to list, but use true controls/state for generating measurements
  Eigen::Matrix4d xn;
  double t = this->ts.back();
  

  this->sys.Step(xn, t, this->xs.back(), u, h);
  this->us.push_back(u);
  this->xs.push_back(xn);
  this->ts.push_back(t + h);

  this->sys.Step(xn, t, this->xos.back(), u, h);
  this->xos.push_back(xn);
  this->uos.push_back(u);
  
  // time-stage 
  int k = this->xs.size()-1;

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
  
  this->Is.resize(this->Is.size()+1);
  assert(k == this->Is.size()-1);

  for (int l = 0; l < this->ls.size(); ++l) {
    const Vector3d &pf = this->ls[l];
    //p.segment<2>(2*extforce + 2*l);
    
    double d = (px - pf).norm();
    if (d < this->dmax) {

      Vector3d z = R.transpose()*(pf - px);
      z /= z.norm();
      z(0) += sqrt(this->cp)*random_normal();
      z(1) += sqrt(this->cp)*random_normal();
      z(2) += sqrt(this->cp)*random_normal();
      z /= z.norm();
      
      // add to feature vector if not observed
      if (!this->observed[l]) {
        if (!this->init) {
          this->p.resize(3);
          this->init= true;
        } else {
          this->p.conservativeResize(this->p.size() + 3);
        }
        const Matrix3d &Rn = this->xs.back().block(0,0,3,3);
        const Vector3d &pxn = this->xs.back().block(0,3,3,1);
         
        // Initialize landmark position from estimated state
        this->p.tail(3) = Rn*this->dmax*z + pxn;
        // Initialize landmark position with its true position
        //p.tail<3>() = pf;
        this->pis.push_back(l);
        this->observed[l] = true;
        this->cis[l] = this->p.size()/3-1;
      }
      
      this->Is[k].push_back(make_pair(l, z)); // add feature l to pose k      
      this->Js[l].push_back(make_pair(k, z)); // add pose k to feature l
      
      //        cout << "k=" << k << " l=" << l << " d=" << d << endl;
    }
  }
}


// Makes the ground truth features to be observed
template <int _nu>
void KinbodyProjTrack<_nu>::MakeTrue()
{
  //  srand (time(NULL));
  //  int nf = (p.size() - extforce*2)/2;
  int i = 0;  
  int bound = 40;
  for (int l =0; l < this->ls.size();) {
    //double a = RND*bound;
    double a = RND*2*M_PI;
    if (a> 0.2 + M_PI && a<0.2 + 1.5*M_PI || a>0.2 && a<0.2 + 0.5*M_PI)
      continue;
    double z = (RND-0.5)*1.5*this->h;
    double r_rand = 3*this->w*(RND-0.5) + this->r; 
    this->ls[l] = Vector3d(r_rand*cos(a), r_rand*sin(a), z);
    //ls[l] = Vector3d(a-bound/2., 5+z, z);
    ++l;
  }
}


// Get the desired state on the track for some time t
template <int _nu>
void KinbodyProjTrack<_nu>::Get(Matrix4d &x, double vd, double t) const
{
  double a = 1.3*2*M_PI*t/this->tf;
  //double a = 20*t/tf;
  x.setIdentity();
  
  x(0,0) = cos(a+M_PI/2.); x(0,1) = -sin(a+M_PI/2.);
  x(1,0) = sin(a+M_PI/2.); x(1,1) = cos(a+M_PI/2.);
 
  x(0,3) = this->r*cos(a); x(1,3) = this->r*sin(a); x(2,3) = 0;
  
  //x(0,3) = 3*a;
}
}

#endif
