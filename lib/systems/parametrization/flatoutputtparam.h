#ifndef GCOP_FLATOUTPUTTPARAM_H
#define GCOP_FLATOUTPUTTPARAM_H

#include "tparam.h"

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
    int np = Dynamic,
    int _ntp = Dynamic> class FlatOutputTparam : public Tparam<T, nx, nu, np, _ntp> {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  

    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, _ntp, 1> Vectorntpd;

	protected:
		/** Evaluate Bezier curve using knots in s using DeCasteljau algorithm
		 * @param u input point where bezier curve is evaluated \in (0,1)
		 * @param start Starting point for recursionP_0
		 * @param end ending point for recursion P_n
		 */
		//VectorXd DeCasteljau(vector<VectorXd> &s, double u, int start, int end);
		VectorXd DeCasteljau(vector<VectorXd> &s, double u, int start, int end)
		{
			if(start == end)
			{
				return s[start];
			}
			else
			{
				return ((1-u)*(this->DeCasteljau(s, u, start, end-1)) + u*(this->DeCasteljau(s,u,start+1,end)));
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
            if(count_knots >= (1 + numberofderivatives) && 
              (!fixfinal || count_knots < numberofknots - (1 + numberofderivatives)))
            {
						  knotsforallderivatives[count_derivatives][count_knots] = s.segment((count_knots-(1+numberofderivatives))*ny,ny);
            }
					}
					else
					{
						knotsforallderivatives[count_derivatives][count_knots] =knotsforallderivatives[count_derivatives-1][count_knots+1] - knotsforallderivatives[count_derivatives-1][count_knots];
					}
					//cout<<"knotsforallderivatives["<<count_derivatives<<"]["<<count_knots<<"]:"<<knotsforallderivatives[count_derivatives][count_knots]<<endl;
					//getchar();//#DEBUG
				}
			}
		}
 
  public:
		/** Constructor
		 * @param sys			System used for manifold
		 * @param numberofknots The number of knots used for bezier curve. This also determines the degree of the curve
		 * @param numberofderivatives_ The number of derivatives needed by the system for evaluating flat outputs
		 */
    FlatOutputTparam(System<T, nx, nu, np> &sys, int ny_, int numberofknots_, int numberofderivatives_ = 0, bool fixfinal_ = false);
    
    bool To(Vectorntpd &s, 
            const vector<double> &ts, 
            const vector<T> &xs, 
            const vector<Vectorcd> &us,
            const Vectormd *p = 0);
    
    bool From(vector<double> &ts, 
              vector<T> &xs, 
              vector<Vectorcd> &us,
              const Vectorntpd &s,
              Vectormd *p = 0);
    
		int numberofderivatives;///< Number of derivatives of flat outputs needed
		int numberofknots;///< Number of knots for bezier curve
		int ny;///< Number of flat outputs
		vector<vector<VectorXd> > knotsforallderivatives;///<Knots for each derivative are computed on the fly
    bool fixfinal;///< should the final flat output be fixed (with necessary flatoutput time derivatives == 0)
  };
  
  template <typename T, int nx, int nu, int np, int _ntp> 
    FlatOutputTparam<T, nx, nu, np, _ntp>::FlatOutputTparam(System<T, nx, nu, np> &sys, int ny_, int numberofknots_, int numberofderivatives_, bool fixfinal_) :  Tparam<T, nx, nu, np, _ntp>(sys, (numberofknots_-(numberofderivatives_+1)-fixfinal_*(numberofderivatives_+1))*ny_),ny(ny_), numberofderivatives(numberofderivatives_), numberofknots(numberofknots_), fixfinal(fixfinal_)  {
			assert(numberofknots > 0);
			assert(ny >0);
			knotsforallderivatives.resize(numberofderivatives+1);
			for(int count = 0;count < numberofderivatives+1; count++)
			{
				int knotsize = numberofknots-count;
				knotsize = knotsize>0?knotsize:0;
				knotsforallderivatives[count].resize(knotsize);
			}
      knotsforallderivatives[0][0] = VectorXd::Zero(ny,1);
  }

  template <typename T, int nx, int nu, int np, int _ntp> 
    bool FlatOutputTparam<T, nx, nu, np, _ntp>::To(Vectorntpd &s, 
                                             const vector<double> &ts, 
                                             const vector<T> &xs, 
                                             const vector<Vectorcd> &us,
																						 const Vectormd *p) {
			//Evaluate the basis matrix at all the points ts:
			int N = us.size();
			double tf = ts.back();
			assert(N+1 > numberofknots);//Assert there are enough points for given number of knots
			MatrixXd basis(ny*(N+1), (numberofknots)*ny);
			basis.setZero();//Initialize to zero
			//Can do inversion separate for each ny. But have to copy back to s separately. #Not very efficient
			VectorXd flatoutputs((N+1)*ny);
			vector<VectorXd> knotvector(numberofknots);
			//Initialize knotvector to zeros to begin with:
			for(int count_knots = 0;count_knots < numberofknots; count_knots++)
			{
				knotvector[count_knots].resize(ny);
				knotvector[count_knots].setZero();
			}


			for(int count_knots =0; count_knots < numberofknots; count_knots++)
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
					basis.block((count_samples)*ny, (count_knots)*ny,ny,ny) = (DeCasteljau(knotvector, (ts[count_samples]/tf), 0,numberofknots-1)).asDiagonal();//The evaluation using P = e_i will be B_n,i(u_i)
				}
			}

			//Evaluate the flat outputs at the given states and controls:
			for(int count_samples = 0; count_samples < N; count_samples++)
			{
				VectorXd flatoutput;
				(this->sys).StateAndControlsToFlat(flatoutput, xs[count_samples], us[count_samples]);
				flatoutputs.segment((count_samples)*ny,ny) = flatoutput;
				cout<<"Flat output["<<count_samples<<"]: "<<flatoutput.transpose()<<endl;
			}
			//Tail or final state:
			{
				VectorXd flatoutput;
				(this->sys).StateAndControlsToFlat(flatoutput, xs[N], us[N-1]);//Copying us[N-1] for us[N] since us[N] does not exist
				flatoutputs.tail(ny) = flatoutput;
				cout<<"Flat output["<<N<<"]: "<<flatoutput.transpose()<<endl;
			}
			//getchar();//#DEBUG
			//cout<<"Basis: "<<endl<<basis<<endl;
			//Set the output knot vector(s) as basis.pseudoinverse()*flatoutputsevaluatedatallsamples
			VectorXd s_temp = ((basis.transpose()*basis).ldlt().solve(basis.transpose()*flatoutputs));
      if(fixfinal)
      { 
        assert(numberofknots-(numberofderivatives+1) >= 0);
        s = s_temp.segment((1+numberofderivatives)*ny, (numberofknots-2*(1+numberofderivatives))*ny);
      }
      else
      {
        s = s_temp.tail((numberofknots-(1+numberofderivatives))*ny);
      }

			cout<<"Error: "<<(basis*s_temp - flatoutputs).squaredNorm()<<endl;

      for(int i = 0; i < 1+numberofderivatives; i++)
      {
        knotsforallderivatives[0][i] = flatoutputs.segment(ny*i, ny);
      }

      if(fixfinal)
      {
        assert(numberofknots-(numberofderivatives+1) >= 0);
        for(int i = 0; i < 1+numberofderivatives; i++)
        {
          knotsforallderivatives[0][numberofknots-(i+1)] = flatoutputs.segment((flatoutputs.size()-ny*(i+1)), ny);
        }
      }
      cout<<"flatoutputs.head(ny) " << flatoutputs.head(ny).transpose() << endl; 
      cout<<"flatoutputs.tail(ny) " << flatoutputs.tail(ny).transpose() << endl; 
      //(this->sys).StateAndControlsToFlat(knotsforallderivatives[0][0], xs[0], us[0]);

			cout<<"s: "<<s.transpose()<<endl;
			cout<<"s_temp: "<<s_temp.transpose()<<endl;
			//s.setZero();
                        return true;
		}
  
  template <typename T, int nx, int nu, int np, int _ntp> 
    bool FlatOutputTparam<T, nx, nu, np, _ntp>::From(vector<double> &ts, 
                                                  vector<T> &xs, 
                                                  vector<Vectorcd> &us,
                                                  const Vectorntpd &s,
																									Vectormd *p) {
			assert(this->ntp == (numberofknots-(numberofderivatives+1)-fixfinal*(numberofderivatives+1))*ny);
			//Evaluate Flat outputs from the knot inputs (s) at all the input times ts
			int N = us.size();
			double tf = ts.back();//Replace this with deltat version #TODO
			vector<VectorXd> flatoutputsandderivatives(numberofderivatives+1);
			Createknotsforallderivatives(s);//Create all the knots before starting evaluation
			for(int count_ts = 0; count_ts <= N; count_ts++)
			{
				double factor_derivative = 1;
				for(int count_derivatives =0; count_derivatives <= numberofderivatives ; count_derivatives++)
				{
					if(numberofknots-count_derivatives>0)//If there are any derivatives to be evaluated then only go further
					{
						flatoutputsandderivatives[count_derivatives] = (factor_derivative)*(DeCasteljau(knotsforallderivatives[count_derivatives], (ts[count_ts]/tf), 0, numberofknots-count_derivatives-1));
				    //cout<<"flatoutputandderivatives["<<count_derivatives<<"]: "<<flatoutputsandderivatives[count_derivatives]<<endl;
					}
					else
					{
						flatoutputsandderivatives[count_derivatives].resize(ny);
						flatoutputsandderivatives[count_derivatives].setZero();
					}
					factor_derivative *= double(numberofknots-count_derivatives)/tf;//tf is here because u = t/tf;
				}
				//getchar();
				//Evaluate system states and controls using the flat outputs and derivatives
        //cout << "flatoutputs::From: " << flatoutputsandderivatives[0].transpose() << endl;
				if(count_ts == 0)
				{
					T xtemp;
					(this->sys).FlatToStateAndControls( xtemp, us[count_ts], flatoutputsandderivatives);
				}
				else if(count_ts == N)
				{
					Vectorcd utemp;//For last point we do not care abt control
					(this->sys).FlatToStateAndControls( xs[count_ts], utemp, flatoutputsandderivatives);
				}
				else
				{
					(this->sys).FlatToStateAndControls( xs[count_ts], us[count_ts], flatoutputsandderivatives);
				}
				//cout<<"us["<<count_ts<<"]: "<<us[count_ts]<<endl;
			}

			//cout<<"s: "<<s.transpose()<<endl;
                        return true;
		}
		
}

#endif
