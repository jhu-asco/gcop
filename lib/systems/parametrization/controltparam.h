#ifndef GCOP_CONTROLTPARAM_H
#define GCOP_CONTROLTPARAM_H

#include "tparam.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  /**
   * General trajectory parametrization functionality. A trajectory can be parametrizes using
   * a set of discrete controls, a continuous control paramertization, a set of discrete states, 
   * a set of discrete flat outputs, etc...
   *
   * The default implementation is the simplest: a set of discrete controls
   *
   * Author: Marin Kobilarov (c) 2005--2013
   */
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic,
    int _ntp = Dynamic> class ControlTparam : public Tparam<T, nx, nu, np, _ntp> {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  

    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, _ntp, 1> Vectorntpd;
    
  public:
    ControlTparam(System<T, nx, nu, np> &sys, const vector<double> &ts);
    
    void To(Vectorntpd &s, 
            const vector<double> &ts, 
            const vector<T> &xs, 
            const vector<Vectorcd> &us,
            const Vectormd *p = 0);
    
    void From(vector<double> &ts, 
              vector<T> &xs, 
              vector<Vectorcd> &us,
              const Vectorntpd &s,
              Vectormd *p = 0);
    
    const vector<double> &tks;  ///< knot times
  };
  
  template <typename T, int nx, int nu, int np, int _ntp> 
    ControlTparam<T, nx, nu, np, _ntp>::ControlTparam(System<T, nx, nu, np> &sys, const vector<double> &tks) :  Tparam<T, nx, nu, np, _ntp>(sys, tks.size()*sys.U.n), tks(tks) {
  }

  template <typename T, int nx, int nu, int np, int _ntp> 
    void ControlTparam<T, nx, nu, np, _ntp>::To(Vectorntpd &s, 
                                             const vector<double> &ts, 
                                             const vector<T> &xs, 
                                             const vector<Vectorcd> &us,
                                             const Vectormd *p) {
    int j =0;//Iterator for tks
    cout<<"s.size: "<<s.size()<<endl;
    for(int i = 0;i < us.size(); ++i)
    {
      if ((ts[i] - tks[j]) >= 0)
      {
        if(j < tks.size())
        {
          s.segment(j*(this->sys.U.n), this->sys.U.n) = us[i];
          cout<<"s["<<j<<"]: "<<us[i].transpose()<<endl;
          j++;
        }
      }
    }
    s.tail(this->sys.U.n) = us[us.size()-1];
    cout<<"s: "<<s.transpose()<<endl;
    //s.setConstant(.001);
  }
  
  template <typename T, int nx, int nu, int np, int _ntp> 
    void ControlTparam<T, nx, nu, np, _ntp>::From(vector<double> &ts, 
                                                  vector<T> &xs, 
                                                  vector<Vectorcd> &us,
                                                  const Vectorntpd &s,
                                                  Vectormd *p) {
    assert(this->ntp == tks.size()*this->sys.U.n);
    
    int ki = 0; // knot index
    for (int i = 0; i < us.size(); ++i) {
      int si = ki*this->sys.U.n;
      assert(si + this->sys.U.n < s.size());
      const Vectorcd &ua = s.segment(si, this->sys.U.n);   // left control
      const Vectorcd &ub = s.segment(si + this->sys.U.n, this->sys.U.n);  // right control

      double &t = ts[i];
      us[i] = ua + (ub - ua)*(t - tks[ki])/(tks[ki+1] - tks[ki]);
      if (t > tks[ki + 1] + 1e-16) {
        ki++;
      }      
    }

    this->sys.reset(xs[0],ts[0]);
    for (int i = 0; i < us.size(); ++i) {
      this->sys.Step_internalinput(xs[i+1], us[i], ts[i+1] - ts[i], p);
      //cout<<"Xs["<<(i+1)<<"]"<<xs[i+1].transpose()<<endl; #DEBUG
      //cout<<"us["<<i<<"]"<<us[i].transpose()<<endl;
    }
  }
}

#endif
