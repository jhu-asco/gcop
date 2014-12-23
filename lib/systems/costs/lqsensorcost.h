#ifndef GCOP_LQSENSORCOST_H
#define GCOP_LQSENSORCOST_H

#include <limits>
#include "lssensorcost.h"
#include <iostream>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  template <typename T = VectorXd, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic,
    int _ng = Dynamic,
    typename Tz = VectorXd,
    int _nz = Dynamic> class LqSensorCost : public LsSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz> {
    public:

  typedef Matrix<double, _ng, 1> Vectorgd;
  typedef Matrix<double, _ng, _nx> Matrixgxd;
  typedef Matrix<double, _ng, _nu> Matrixgud;
  typedef Matrix<double, _ng, _np> Matrixgpd;
  typedef Matrix<double, _ng, _nz> Matrixgzd;

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;

  typedef Matrix<double, _nx, _nx> Matrixnd;
  typedef Matrix<double, _nx, _nu> Matrixncd;
  typedef Matrix<double, _nu, _nx> Matrixcnd;
  typedef Matrix<double, _nu, _nu> Matrixcd;
  
  //  typedef Matrix<double, Dynamic, 1> Vectormd;
  typedef Matrix<double, _np, _np> Matrixmd;
  typedef Matrix<double, _nx, _np> Matrixnmd;
  typedef Matrix<double, _np, _nx> Matrixmnd;
      
  typedef Matrix<double, _nz, 1> Vectorrd;
  typedef Matrix<double, _nz, _nz> Matrixrd;
  typedef Matrix<double, _nz, _nx> Matrixrnd;
  typedef Matrix<double, _nz, _nu> Matrixrcd;
  typedef Matrix<double, _nz, _np> Matrixrmd;

    /**
     * Linear-quadratic cost on a sensor manifold. Use this constructor for dynamic-size
     * control problem, i.e. LqCost<T>(X, U, tf, xf, ...)
     * @param sys provides system manifold
     * @param sensor provides sensor manifold
     * @param tf final time
     * @param diag whether the cost matrices are diagonal?
     */
    LqSensorCost(System<T, _nx, _nu, _np> & sys, Manifold<Tz, _nz>  &Z, bool diag = true);

    /**
     * Updates the square root gains, i.e. sets Qsqrt = Q.sqrt(), etc...
     *
     */
    void UpdateGains(); 

    double L(double t, const Tz &z, const Vectornd &w,
                   const Vectormd &p, double h, int sensor_index, 
                   Vectornd *Lw = 0, Matrixnd* Lww = 0,
                   Vectormd *Lp = 0, Matrixmd* Lpp = 0,
                   Vectorrd *Lz = 0, Matrixrd* Lzz = 0,
                   Matrixmnd *Lpw = 0, Matrixrnd* Lzw = 0, Matrixrmd* Lzp = 0);  

    double Lp(const Vectormd &p,
                   Vectormd *Lp = 0, Matrixmd *Lpp = 0);

    bool Res(Vectorgd &g, 
                     double t, const Tz &z,
                     const Vectornd &w, const Vectormd &p, 
                     double h, int sensor_index,
                     Matrixgxd *dgdw = 0, Matrixgpd *dgdp = 0,
                     Matrixgzd *dgdz = 0);    

    bool Resp(Vectormd &gp, 
                   const Vectormd &p,
                   Matrixmd *dgdp = 0);    
    /**
     * Set observed measurements for the cost
     * @param zs observed measurements
     * @param mup mean parameters (Prior for parameters)
     */
    void SetReference(const vector<Tz> *zs, Vectormd *mup = 0);        
    
    //const Tz &zf; ///< reference to a desired final state
    
    bool diag;   ///< are the Q, Qf, and R matrices diagonal (true by default)
    
    Matrixrd R;       ///< inverse of measurement covariance matrix
    Matrixnd S;       ///< inverse of noise covariance matrix
    Matrixmd P;       ///< inverse of parameter covariance matrix (set to zero if nothing is known about parameters before hand)

    Matrixnd Ssqrt;       ///< S sqrt
    Matrixrd Rsqrt;       ///< R sqrt
    Matrixmd Psqrt;       ///< P sqrt

    const vector<Tz> *zs;         ///< measured sensor data
    Vectormd *mup;          ///< Prior for parameter data

    protected:

    Vectorrd dz;      ///< measurement error (as a tangent vector)
    Vectormd dp;      ///< parameter error
    };
  
    template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz>
    LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::LqSensorCost(
                                                    System<T, _nx, _nu, _np> &sys, 
                                                    Manifold<Tz, _nz>  &Z,
                                                    bool diag) : 
    LsSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>(sys, Z, sys.X.n + Z.n + sys.P.n), diag(diag) {
    
    if (_nz == Dynamic || _nx == Dynamic || _np == Dynamic) {
      R.resize(Z.n, Z.n);
      S.resize(sys.X.n, sys.X.n);
      P.resize(sys.P.n, sys.P.n);

      Rsqrt.resize(Z.n, Z.n);
      Ssqrt.resize(sys.X.n, sys.X.n);
      Psqrt.resize(sys.P.n, sys.P.n);

      dp.resize(sys.P.n);
    }

    R.setIdentity();
    S.setIdentity();
    P.setZero();//Initial guess for parameters should be no prior

    UpdateGains();

    zs = 0;
    mup = 0;
  }
  
  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz>
    void LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::UpdateGains() {
    Rsqrt = R.array().sqrt();
    Ssqrt = S.array().sqrt();
    Psqrt = P.array().sqrt();
  }


  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz>
    void LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::SetReference(const vector<Tz> *zs, Vectormd *mup) {
    this->zs = zs;
    //this->mup = mup;
    if(mup)
    {
      this->mup = new VectorXd(mup->size());
      *(this->mup) = *mup;//Copy the vector over
    }
    //Add assert statments to make sure there are enough measurements and parameters#TODO
  }

  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz>
    bool LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::Res(Vectorgd &g, 
                     double t, const Tz &z,
                     const Vectornd &w, const Vectormd &p, 
                     double h,int sensor_index,
                     Matrixgxd *dgdw, Matrixgpd *dgdp,
                     Matrixgzd *dgdz)
    {
      assert(g.size() == this->sys.X.n + this->sys.P.n + this->Z.n);

      int &nw = this->sys.X.n;
      int &nz = this->Z.n;
      int &np = this->sys.P.n;

      assert(h > 0);
      //int k = round(t/h);
      if (zs) {
        assert(sensor_index < zs->size());
        this->Z.Lift(dz, (*zs)[sensor_index], z); // difference (on a vector space we have dx = x - xf)
        //std::cout<<"dz: "<<dz.transpose()<<"\t"<<t<<"\t"<<h<<endl;
        //std::cout<<"zs["<<k<<"]: "<<(*zs)[k].transpose()<<endl;
        //std::cout<<"z["<<k<<"]: "<<z.transpose()<<endl;
      } 
      //cout<<"zs["<<sensor_index<<"]: "<<(*zs)[sensor_index].transpose()<<"\t"<<z.transpose()<<endl;
      //getchar();
      assert(!std::isnan(dz[0]));

      if (diag)
      {
        g.head(nw) = Ssqrt.diagonal().cwiseProduct(sqrt(h)*w);
        g.segment(nw, nz) = Rsqrt.diagonal().cwiseProduct(sqrt(h)*dz);
  //     cout<<"w: "<<w.transpose()<<endl;
      }
      else
      {
        g.head(nw) = Ssqrt*(sqrt(h)*w);
        g.segment(nw, nz) = Rsqrt*(sqrt(h)*dz);
      }
      g.tail(np).setZero();
      //cout<<"g: "<<g.transpose()<<endl;

      return true;
    }

  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz>
    bool LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::Resp(Vectormd &gp, 
                                                            const Vectormd &p,
                                                            Matrixmd *dgdp)
    {
      if (mup) {
        dp = p - *mup;//Error
      }
      else
        dp  = p;

      if(diag)
      {
        gp = Psqrt.diagonal().cwiseProduct(dp);
      }
      else
      {
        gp = Psqrt*dp;
      }
    }

  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz>
    double LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::L(double t, const Tz &z, const Vectornd &w,
                                                           const Vectormd &p, double h, int sensor_index,
                                                           Vectornd *Lw, Matrixnd* Lww,
                                                           Vectormd *Lp, Matrixmd* Lpp,
                                                           Vectorrd *Lz, Matrixrd* Lzz,
                                                           Matrixmnd *Lpw, Matrixrnd* Lzw, Matrixrmd* Lzp)
    {
      int k = (int)(t/h);

      if (zs) {
        assert(sensor_index < zs->size());
        this->Z.Lift(dz, (*zs)[sensor_index], z); // difference (on a vector space we have dx = x - xf)
      } 
      assert(!std::isnan(dz[0]));

      if (Lz)
        if (diag)
          *Lz = R.diagonal().cwiseProduct(h*dz);
        else
          *Lz = R*(h*dz);

      if (Lzz) 
        *Lzz = h*R;

      if (Lw)
        if (diag)
          *Lw = S.diagonal().cwiseProduct(h*w);
        else
          *Lw = S*(h*w);        

      if (Lww)
        *Lww = h*S;

      if(Lp)
        Lp->setZero();

      if(Lpp)
        Lpp->setZero();

      if (Lzw)
        Lzw->setZero();

      if(Lzp)
        Lzp->setZero();

      if(Lpw)
        Lpw->setZero();

      if (diag)
        return h*(dz.dot(R.diagonal().cwiseProduct(dz))/2 + w.dot(S.diagonal().cwiseProduct(w)))/2;
      else
        return h*(dz.dot(R*dz) + w.dot(S*w))/2;
      return 0;
    }  

  template <typename T, int _nx, int _nu, int _np, int _ng, typename Tz, int _nz>
    double LqSensorCost<T, _nx, _nu, _np, _ng, Tz, _nz>::Lp(const Vectormd &p,
                                                            Vectormd *Lp, Matrixmd *Lpp)
    {
      if(mup)
        dp = p - *mup;
      else
        dp = p;

      if (Lp)
        if (diag)
          *Lp = P.diagonal().cwiseProduct(dp);
        else
          *Lp = P*dp;        

      if (Lpp)
        *Lpp = P;

      if (diag)
        return dp.dot(P.diagonal().cwiseProduct(dp))/2;
      else
        return dp.dot(P*dp)/2;
    }

}


#endif
