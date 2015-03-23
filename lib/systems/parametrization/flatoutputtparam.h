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
	 * TODO If ny Dynamic somehow resize all the knotvector sizes
	 *
   * Author: Marin Kobilarov (c) 2005--2013
   * Author2: Gowtham Garimella
   */
  template <typename T, 
    int nx, 
    int nu,
		int ny,
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
		 * @param numberofderivatives_ The number of derivatives needed by the system for evaluating flat outputs
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
			//Evaluate the basis matrix at all the points ts:
			int N = us.size();
			double tf = ts.back();
			MatrixXd basis(ny*(N+1), numberofknots);
			//Can do inversion separate for each ny. But have to copy back to s separately. #Not very efficient
			VectorXd flatoutputs((N+1)*ny);
			vector<Vectoryd> knotvector(knotsize);
			//Set knotvector to zeros to begin with:
			for(int count_knots = 0;count_knots < knotsize; count_knots++)
				knotvector[count_knots].setZero();

			for(int count_samples =0; count_samples < N; count_samples++)
			{
				//Evaluate the flat outputs at the given states and controls:
				(this->sys).StateAndControlsToFlat(flatoutputs.segment<ny,1>(count_samples*ny), xs[count_samples], us[count_samples]);
			}
			//Tail or final state:
			(this->sys).StateAndControlsToFlat(flatoutputs.tail<ny,1>(), xs[N], us[N-1]);//Copying us[N-1] for us[N] since us[N] does not exist

			for(int count_knots =0; count_knots < knotsize; count_knots++)
			{
				//Set the Knots to be e_i and evaluate across all the samples to find the basis matrix
				if(count_knots > 0)
				{
					knotvector[count_knots-1].setZero();
				}
				knotvector[count_knots].setConstant(1.0);

				for(int count_samples =0; count_samples <= N; count_samples++)
				{
					//Evaluate the Bezier curve at given ts using the above knot vector:
					basis.block<ny,1>(count_samples*ny, count_knots) = DeCasteljau(knotvector, (ts[count_samples]/tf), 0,numberofknots-1);//The evaluation using P = e_i will be B_n,i(u_i)
				}
			}
			//Set the output knot vector(s) as basis.pseudoinverse()*flatoutputsevaluatedatallsamples
			s = (basis.transpose()*basis).ldlt().solve(basis.transpose()*flatoutputs);
			cout<<"Error: "<<(basis*s - flatoutputs).squaredNorm()<<endl;

			cout<<"s: "<<s.transpose()<<endl;
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
			double tf = ts.back();
			vector<Vectoryd> flatoutputsandderivatives(numberofderivatives+1);
			Createknotsforallderivatives(s);//Create all the knots before starting evaluation
			for(int count_ts = 0; count_ts <= N; count_ts++)
			{
				for(int count_derivatives =0; count_derivatives <= numberofderivatives ; count_derivatives++)
				{
					if(numberofknots-count_derivatives>0)//If there are any derivatives to be evaluated then only go further
					{
						flatoutputsandderivatives(count_derivatives) = DeCasteljau(knotsforallderivatives(count_derivatives), (ts[count_ts]/tf), 0, numberofknots-1);
					}
					else
					{
						flatoutputsandderivatives(count_derivatives).setZero();
					}
				}
				//Evaluate system states and controls using the flat outputs and derivatives
				if(count_ts < N)
					(this->sys).FlatToStateAndControls( xs[count_ts], u[count_ts], flatoutputsandderivatives);
				else
				{
					Vectorcd utemp;//For last point we do not care abt control
					(this->sys).FlatToStateAndControls( xs[count_ts], utemp, flatoutputsandderivatives);
				}
			}
		}
}

#endif
