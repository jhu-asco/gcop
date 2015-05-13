#ifndef GCOP_SPLINETPARAM_H
#define GCOP_SPLINETPARAM_H

#include "tparam.h"
#include <unsupported/Eigen/Splines>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  /**
   * This version implements the non uniform spline interpolation of the data. It also implements least square fitting of the given control data using spline
   *
   * Author: Marin Kobilarov (c) 2005--2013
   * Author2: Gowtham Garimella
   */
  template <typename T, 
    int nx, 
    int nu,
    int np = Dynamic,
    int _ntp = Dynamic> class SplineTparam : public Tparam<T, nx, nu, np, _ntp> {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  

    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, _ntp, 1> Vectorntpd;

    typedef Spline<double, nu> ControlSpline;
 
  public:
		/** Constructor
		 * @param sys			System used for manifold
		 * @param tks			Knot times for spline
		 * @param degree	Spline degree
		 */
    SplineTparam(System<T, nx, nu, np> &sys, const VectorXd tks, int degree = 2);//By default use quadratic spline
    
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
    
    VectorXd tks;  ///< control times
    double tf;///< Final time for spline
    int degree; ///< Degree of the spline p = (m - n - 1) where m is knot vector size(used by bspline to find the basis etc); n is the control vector size (tks size)
  };
  
  template <typename T, int nx, int nu, int np, int _ntp> 
    SplineTparam<T, nx, nu, np, _ntp>::SplineTparam(System<T, nx, nu, np> &sys, const VectorXd tks, int degree) :  Tparam<T, nx, nu, np, _ntp>(sys, tks.size()*sys.U.n), tks(tks), degree(degree) {
      //normalize tks from 0 to 1
      assert(tks.size() > degree+1);//Make sure there are enough points
      tf = this->tks(tks.size()-1);//Last element
      //Assuming  increasing order normalize so that last element is 1.0:
      (this->tks) /= tf;
      this->tks(tks.size()-1) = 1.0;//Redundancy
      cout<<"tks: "<<tks.transpose()<<endl;
  }

  template <typename T, int nx, int nu, int np, int _ntp> 
    void SplineTparam<T, nx, nu, np, _ntp>::To(Vectorntpd &s, 
                                             const vector<double> &ts, 
                                             const vector<T> &xs, 
                                             const vector<Vectorcd> &us,
                                             const Vectormd *p) {
      //Find control points(uks) given us
      int m = degree + 1 + tks.size();//Number of knots
      int nofsegments = m - 2*(degree+1)+1;
      assert(nofsegments > 0);
      //Create knot vector:
      VectorXd knotVector(m);
      KnotAveraging(tks, degree, knotVector);
      //cout<<"knotvector: "<<knotVector.transpose()<<endl;
      //cout<<"tks: "<<tks.transpose()<<endl;

      //Initialize variables for least square fit; Fitting each control dimension separately
      MatrixXd A(us.size(), tks.size());//Basis Matrix
      MatrixXd Asquare(tks.size(), tks.size());//Basis Matrix
      VectorXd c(us.size());//RHS 
      A.setZero();//Initialize A
      //VectorXd b(tks.size());//control points

      //Find Basis Matrix to find the control points:
      double tdiff = ts.back() - ts.front();
      for(int ind = 0;ind < us.size(); ind++)
      {
        double ts_normalized = (ts[ind] - ts[0])/tdiff;
        int index = ControlSpline::Span(ts_normalized, degree, knotVector);
        VectorXd basisFcns = ControlSpline::BasisFunctions(ts_normalized, degree, knotVector);
       // cout<<"Index: "<<index<<endl;
        //cout<<"BasisFcns: "<<basisFcns.transpose()<<endl;
        A.block(ind, index-degree, 1, degree+1) = basisFcns.transpose();//Create Basis Matrix
      }
      Asquare = (A.transpose()*A);
      //cout<<"A: "<<endl<<A<<endl;
      //cout<<"Asquare: "<<endl<<Asquare<<endl;
      //getchar();
      for(int ind = 0; ind < nu; ind++)
      {
        for(int uind = 0; uind < us.size(); uind++)
        {
          c(uind) = us[uind](ind);
        }
        VectorXd b = Asquare.ldlt().solve(A.transpose()*c);
        cout<<"Error: "<<(A*b - c).squaredNorm()<<endl;
      //  cout<<"c: "<<c.transpose()<<endl;
        //cout<<"b: "<<b.transpose()<<endl;
        //Copy the elements back into vector s:
        for(int sind = 0; sind < tks.size(); sind++)
        {
          s(sind*nu + ind) = b(sind);
        }
      }
      cout<<"s: "<<s.transpose()<<endl;
      //Verify if this s is good:
      /*{
        MatrixXd controlMatrix(this->sys.U.n, tks.size());
        for(int cs = 0; cs < tks.size(); cs++)
        {
          controlMatrix.col(cs) = s.segment(cs*nu, nu);//Fill Control Matrix
        }
        ControlSpline cspline = SplineFitting<ControlSpline>::Interpolate(controlMatrix, degree, tks);
        for (int i = 0; i < us.size(); ++i) {
          Vectorcd usbar = cspline((ts[i] - ts[0])/tdiff);
          cout<<"us_pred["<<i<<"]: "<<usbar.transpose()<<"ts: "<<ts[i]<<endl;
          cout<<"us_actual["<<i<<"]: "<<us[i].transpose()<<"ts: "<<ts[i]<<endl;
        }
        getchar();
      }
      */

      //s.setZero();
    }
  
  template <typename T, int nx, int nu, int np, int _ntp> 
    void SplineTparam<T, nx, nu, np, _ntp>::From(vector<double> &ts, 
                                                  vector<T> &xs, 
                                                  vector<Vectorcd> &us,
                                                  const Vectorntpd &s,
                                                  Vectormd *p) {
 

    assert(this->ntp == tks.size()*this->sys.U.n);

    //Create knot matrix:
    MatrixXd controlMatrix(this->sys.U.n, tks.size());
    for(int cs = 0; cs < tks.size(); cs++)
    {
      controlMatrix.col(cs) = s.segment(cs*nu, nu);//Fill Control Matrix
    }
    //cout<<"s: "<<s.transpose()<<endl;//Input control points
    //cout<<"Control Matrix: "<<endl<<controlMatrix<<endl;
    ControlSpline cspline = SplineFitting<ControlSpline>::Interpolate(controlMatrix, degree, tks);
    double tdiff = ts.back() - ts.front();
    
    this->sys.Reset(xs[0],ts[0]);
    for (int i = 0; i < us.size(); ++i) {
      us[i] = cspline((ts[i] - ts[0])/tdiff);
      //cout<<"ts["<<i<<"]: "<<ts[i]<<endl;
      //cout<<"us["<<i<<"]: "<<us[i].transpose()<<"ts: "<<ts[i]<<endl;
      this->sys.Step(xs[i+1], us[i], ts[i+1] - ts[i], p);
      //cout<<"Xs["<<(i+1)<<"]"<<xs[i+1].transpose()<<endl; #DEBUG
      //cout<<"us["<<i<<"]"<<us[i].transpose()<<endl;
    }
    //getchar();
  }
}

#endif
