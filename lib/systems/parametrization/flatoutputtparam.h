#ifndef GCOP_FLATOUTPUTTPARAM_H
#define GCOP_FLATOUTPUTTPARAM_H

#include "tparam.h"
#include <unsupported/Eigen/Splines>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  /**
   *
	 * This is a flat output parametrization class. The template vector ny gives the size of flat output
	 *
   * Author: Marin Kobilarov (c) 2005--2013
   * Author2: Gowtham Garimella
   */
  template <typename T, 
    int nx, 
    int nu,
		int ny = Dynamic,
    int np = Dynamic,
    int _ntp = Dynamic> class FlatOutputTparam : public Tparam<T, nx, nu, np, _ntp> {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, ny, 1> Vectoryd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  

    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, _ntp, 1> Vectorntpd;

    typedef Spline<double, nu> ControlSpline;

	protected:
		/** Evaluate Bezier curve using knots in s using DeCasteljau algorithm
		 * @param u input point where bezier curve is evaluated \in (0,1)
		 * @param start Starting point for recursionP_0
		 * @param end ending point for recursion P_n
		 */
		Vectoryd DeCasteljau(vector<Vectoryd> &s, double u, int start, int end)
		{
			if(start == end)
			{
				return s(start);
			}
			else
			{
				return ((1-u)*DeCasteljau(s, u, start, end-1) + u*Decasteljau(s,u,start+1,end));
			}
		}

		/** This function evaluates all the knots Dk_i for all the derivatives needed
		 * This is dont recursively as noted in  	http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html
		 * @param s The input knots
		 */
		void Createknotsforallderivatives(const Vectorntpd &s)
		{
			for(int count_derivatives = 0; count_derivatives <= numberofderivatives; count_derivatives++)
			{
				int knotsize = (numberofknots-count_derivatives)>0?(numberofknots-count_derivatives):0;
				for(int count_knots = 0;count_knots < knotsize; count_knots++)
				{
					if(count_derivatives == 0)
					{
						knotsforallderivatives[count_derivatives][count_knots] = s.segment<ny>(count_knots*ny);
					}
					else
					{
						knotsforallderivatives[count_derivatives][count_knots] =knotsforallderivatives[count_derivatives-1][count_knots+1] - knotsforallderivatives[count_derivatives-1][count_knots];
					}
				}
			}
		}
 
  public:
		/** Constructor
		 * @param sys			System used for manifold
		 * @param numberofknots The number of knots used for bezier curve. This also determines the degree of the curve
		 */
    FlatOutputTparam(System<T, nx, nu, np> &sys, int numberofknots_, int numberofderivatives_ = 0);
    
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
    
    //VectorXd tks;  ///< control times
		int numberofderivatives;///< Number of derivatives of flat outputs needed
		int numberofknots;///< Number of knots for bezier curve
		vector<vector<Vectoryd> > knotsforallderivatives;///<Knots for each derivative are computed on the fly
  };
  
  template <typename T, int nx, int nu, int np, int _ntp> 
    FlatOutputTparam<T, nx, nu, ny, np, _ntp>::FlatOutputTparam(System<T, nx, nu, np> &sys, int numberofknots_, int numberofderivatives_) :  Tparam<T, nx, nu, np, _ntp>(sys, numberofknots_*ny), tks(tks), degree(degree), numberofderivatives(numberofderivatives_), numberofknots(numberofknots_)  {
			assert(numberofknots > 0);
			knotsforallderivatives.resize(numberofderivatives+1);
			for(int count = 0;count < numberofderivatives+1; count++)
			{
				int knotsize = numberofknots-count;
				knotsize = knotsize>0?knotsize:0;
				knotsforallderivatives.resize(knotsize);
			}
  }

  template <typename T, int nx, int nu, int np, int _ntp> 
    void FlatOutputTparam<T, nx, nu, np, _ntp>::To(Vectorntpd &s, 
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
    void FlatOutputTparam<T, nx, nu, np, _ntp>::From(vector<double> &ts, 
                                                  vector<T> &xs, 
                                                  vector<Vectorcd> &us,
                                                  const Vectorntpd &s,
                                                  Vectormd *p) {
 
    assert(this->ntp == numberofknots*ny);
		//Evaluate Flat outputs from the knot inputs (s) at all the input times ts
		int N = us.size();
		vector<Vectoryd> flatoutputsandderivatives(numberofderivatives+1);
		Createknotsforallderivatives(s);//Create all the knots before starting evaluation
		for(int count_ts = 0; count_ts <= N; count_ts++)
		{
			for(int count_derivatives =0; count_derivatives <= numberofderivatives ; count_derivatives++)
			{
				if(numberofknots-count_derivatives>0)
				{
					flatoutputsandderivatives(count_derivatives) = DeCasteljau(knotsforallderivatives(count_derivatives), 0, numberofknots-1);
				}
				else
				{
					flatoutputsandderivatives(count_derivatives).setZero();
				}
			}
			//Evaluate system states and controls using the flat outputs and derivatives
		}


    for(int cs = 0; cs < tks.size(); cs++)
    {
      controlMatrix.col(cs) = s.segment(cs*nu, nu);//Fill Control Matrix
    }
    //cout<<"s: "<<s.transpose()<<endl;//Input control points
    //cout<<"Control Matrix: "<<endl<<controlMatrix<<endl;
    ControlSpline cspline = SplineFitting<ControlSpline>::Interpolate(controlMatrix, degree, tks);
    double tdiff = ts.back() - ts.front();
    
    this->sys.reset(xs[0],ts[0]);
    for (int i = 0; i < us.size(); ++i) {
      us[i] = cspline((ts[i] - ts[0])/tdiff);
      //cout<<"ts["<<i<<"]: "<<ts[i]<<endl;
      //cout<<"us["<<i<<"]: "<<us[i].transpose()<<"ts: "<<ts[i]<<endl;
      this->sys.Step_internalinput(xs[i+1], us[i], ts[i+1] - ts[i], p);
      //cout<<"Xs["<<(i+1)<<"]"<<xs[i+1].transpose()<<endl; #DEBUG
      //cout<<"us["<<i<<"]"<<us[i].transpose()<<endl;
    }
    //getchar();
  }
}

#endif
