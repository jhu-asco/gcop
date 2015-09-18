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
		/** Perform (N choose K) operation
		 * @param n N
		 * @param k K
		 */
    double NChooseK(const int n, const int k)
    {
      if(k>n/2)
      {
        return NChooseK(n,n-k);
      }
      else if(k==1)
      {
        return n;
      }
      else
      {
        double c = 1;
        for(int i = 1;i<=k;i++)
        {
          c *= (((double)n-k+i)/((double)i));
        }
        return std::round(c);
      }
    }
		/** Evaluate Bezier curve using knots in s using binomial coefficients algorithm
		 * @param u input point where bezier curve is evaluated \in (0,1)
		 * @param start Starting point for P_0
		 * @param end ending point for P_n
		 */
    VectorXd BinomialInterp(vector<VectorXd> &s, double u, int start, int end)
    { 
      assert(s.size() > 0);
      VectorXd ret(s[0].size());
      ret.setZero();
      for(int i = 0; i <= end; i++)
      {
        ret += NChooseK(end,i)*pow(1-u, end-i)*pow(u,i)*s.at(i);
      }
      return ret;
    }
		/** Evaluate Bezier curve using knots in s using DeCasteljau algorithm
		 * @param u input point where bezier curve is evaluated \in (0,1)
		 * @param start Starting point for recursionP_0
		 * @param end ending point for recursion P_n
		 */
		VectorXd DeCasteljau(vector<VectorXd> &s, double u, int start, int end)
		{
			//if(start == end)
			if(end == 0)
			{
				return s[start];
			}
			else
			{
				return ((1-u)*(this->DeCasteljau(s, u, start, end-1)) + u*(this->DeCasteljau(s,u,start+1,end-1)));
			}
		}

    Eigen::MatrixXd CreateDerivativeBasis(double u, double tf, int deriv_num)
    {
			MatrixXd basis(ny, (numberofknots)*ny);
			basis.setZero();//Initialize to zero
			vector<VectorXd> knotvector(numberofknots);
			for(int count_knots = 0;count_knots < numberofknots; count_knots++)
			{
				knotvector[count_knots].resize(ny);
				knotvector[count_knots].setZero();
			}

			for(int count_knots =0; count_knots < numberofknots; count_knots++)
			{
        vector<vector<VectorXd>> knot_derivatives;
        if(count_knots > 0)
        {
				  knotvector[count_knots-1].setZero();
        }
        knotvector[count_knots].setConstant(1.0);
        
        knot_derivatives.resize(numberofderivatives+1);
        for(int count = 0;count < numberofderivatives+1; count++)
        {
          int knotsize = numberofknots-count;
          knotsize = knotsize>0?knotsize:0;
          knot_derivatives[count].resize(knotsize);
        }
        knot_derivatives[0] = knotvector;
        CreateKnotDerivativesFromKnots(knot_derivatives); 

        double factor_derivative = 1.0;
        for(int count_derivs =0; count_derivs < deriv_num; count_derivs++)
        {
					factor_derivative *= double(numberofknots-count_derivs)/tf;//tf is here because u = t/tf;
        }
				//for(int count_derivs =0; count_derivs < fixed_derivatives+1; count_derivs++)
				//{
					//Evaluate the Bezier curve at given ts using the above knot vector:
				  basis.block(0, count_knots*ny, ny, ny) = 
              (factor_derivative)*(BinomialInterp(knot_derivatives[deriv_num], 
              u, 0, numberofknots-deriv_num-1)).asDiagonal();
				//}
			}
      return basis;
    }

    Eigen::MatrixXd CreateAllDerivativeBasis(double u, double tf)
    {
			MatrixXd basis(ny*(fixed_derivatives+1), (numberofknots)*ny);
			basis.setZero();//Initialize to zero
			vector<VectorXd> knotvector(numberofknots);
			for(int count_knots = 0;count_knots < numberofknots; count_knots++)
			{
				knotvector[count_knots].resize(ny);
				knotvector[count_knots].setZero();
			}

			for(int count_knots =0; count_knots < numberofknots; count_knots++)
			{
        vector<vector<VectorXd>> knot_derivatives;
        if(count_knots > 0)
        {
				  knotvector[count_knots-1].setZero();
        }
        knotvector[count_knots].setConstant(1.0);
        
        knot_derivatives.resize(numberofderivatives+1);
        for(int count = 0;count < numberofderivatives+1; count++)
        {
          int knotsize = numberofknots-count;
          knotsize = knotsize>0?knotsize:0;
          knot_derivatives[count].resize(knotsize);
        }
        knot_derivatives[0] = knotvector;
        CreateKnotDerivativesFromKnots(knot_derivatives); 

        double factor_derivative = 1.0;
				for(int count_derivs =0; count_derivs < fixed_derivatives+1; count_derivs++)
				{
					//Evaluate the Bezier curve at given ts using the above knot vector:
				  basis.block((count_derivs)*ny, (count_knots)*ny,ny,ny) = 
              (factor_derivative)*(BinomialInterp(knot_derivatives[count_derivs], 
              u, 0, numberofknots-count_derivs-1)).asDiagonal();
					factor_derivative *= double(numberofknots-count_derivs)/tf;//tf is here because u = t/tf;
				}
			}
      return basis;
    }

    void CreateKnotDerivativesFromKnots(vector<vector<VectorXd> >& knot_derivatives)
    {
			for(int count_derivatives = 1; count_derivatives <= numberofderivatives; count_derivatives++)
			{
				int knotsize = (numberofknots-count_derivatives)>0?(numberofknots-count_derivatives):0;
				for(int count_knots = 0;count_knots < knotsize; count_knots++)
				{
				knot_derivatives[count_derivatives][count_knots] =
          knot_derivatives[count_derivatives-1][count_knots+1] 
          - knot_derivatives[count_derivatives-1][count_knots];
				}
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
					  knotsforallderivatives[count_derivatives][count_knots] = 
              s.segment((count_knots)*ny,ny);
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

    Eigen::MatrixXd CalcNullspace(const Eigen::MatrixXd& m)
    {
      Eigen::JacobiSVD<MatrixXd> svd(m ,Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::MatrixXd nullspace(m.cols(), svd.matrixV().cols()-svd.nonzeroSingularValues());
      nullspace = svd.matrixV().rightCols(svd.matrixV().cols()-svd.nonzeroSingularValues());
      return nullspace;
    }
 
    Eigen::MatrixXd pinv(const Eigen::MatrixXd& m)
    {
      const double epsilon = std::numeric_limits<double>::epsilon();//1e-6;
      Eigen::JacobiSVD<MatrixXd> svd(m ,Eigen::ComputeThinU | Eigen::ComputeThinV);
	    double tolerance = 
        epsilon * std::max(m.cols(), m.rows()) *svd.singularValues().array().abs()(0);
	    return svd.matrixV() *  
        (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
    }
  public:
		/** Constructor
		 * @param sys			System used for manifold
		 * @param numberofknots The number of knots used for bezier curve. This also determines the degree of the curve
		 * @param numberofderivatives_ The number of derivatives needed by the system for evaluating flat outputs
		 */
    FlatOutputTparam(System<T, nx, nu, np> &sys, int ny_, int numberofknots_, int numberofderivatives_ = 0, int fixed_derivatives_ = 0, bool fixfinal_ = false);
    
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
    bool fixfinal;///< should the final flat output be fixed (with necessary flatoutput time derivatives also fixed)
    MatrixXd fixed_condition_nullspace;
    MatrixXd fixed_condition_basis;
    VectorXd fixed_condition_knot_vec;
    int fixed_derivatives;
    vector<vector<VectorXd> > flatoutputsandderivatives_all;
    Eigen::MatrixXd W;
  };
  
  template <typename T, int nx, int nu, int np, int _ntp> 
    FlatOutputTparam<T, nx, nu, np, _ntp>::FlatOutputTparam(System<T, nx, nu, np> &sys, int ny_, 
      int numberofknots_, int numberofderivatives_, int fixed_derivatives_, bool fixfinal_) :  
      Tparam<T, nx, nu, np, _ntp>(sys, 
        (numberofknots_-(fixfinal_?2:1)*(1+fixed_derivatives_))*ny_),
        ny(ny_), numberofderivatives(numberofderivatives_), numberofknots(numberofknots_), 
        fixfinal(fixfinal_)  {
			assert(numberofknots >= (fixfinal?2:1)*(1+fixed_derivatives_));
			assert(ny >0);
      assert(numberofderivatives_ >= fixed_derivatives_);
      fixed_derivatives = fixed_derivatives_; 
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
			MatrixXd basis(2*ny*(N+1), (numberofknots)*ny);
			basis.setZero();//Initialize to zero
			//Can do inversion separate for each ny. But have to copy back to s separately. #Not very efficient
			VectorXd flatoutputs(2*(N+1)*ny);
      flatoutputs.setZero();
			vector<VectorXd> knotvector(numberofknots);

      // Numerically differeniate first derivative of controls in case it is needed
      vector<vector<Vectorcd> > us_and_ders;
      us_and_ders.resize(N);
      for(int u_count = 0; u_count < us.size(); u_count++)
      {
        us_and_ders[u_count].push_back(us[u_count]);
        if(u_count >= 1 && u_count < us.size()-1)
        {
          us_and_ders[u_count].push_back((us[u_count+1]-us[u_count])/(ts[u_count+1]-ts[u_count]));
        }
        else
        {
          Vectorcd u_zero;
          u_zero.setZero();
          us_and_ders[u_count].push_back(u_zero);
        }
      }

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
					basis.block((count_samples)*ny, (count_knots)*ny,ny,ny) = (BinomialInterp(knotvector, (ts[count_samples]/tf), 0,numberofknots-1)).asDiagonal();//The evaluation using P = e_i will be B_n,i(u_i)
				}
			}
      
      
      // Do minimum snap!      
      for(int count_samples = 0; count_samples <= N; count_samples++)
      {
        basis.block(ny*(N+1 + count_samples), 0, ny, (numberofknots)*ny) = 
          CreateDerivativeBasis(ts[count_samples]/tf, tf, 4);
      }
      flatoutputs.tail(ny*(N+1)).setZero();
      
      /*
      // Fit acceleration
      for(int count_samples = 0; count_samples <= N; count_samples++)
      {
        basis.block(ny*(N+1 + count_samples), 0, ny, (numberofknots)*ny) = 
          CreateDerivativeBasis(ts[count_samples]/tf, tf, 2);
      }
      */
			//Evaluate the flat outputs at the given states and controls:
			for(int count_samples = 0; count_samples < N; count_samples++)
			{
				vector<VectorXd> flatoutput;
				(this->sys).StateAndControlsToFlatAndDerivatives(flatoutput, xs[count_samples], 
          us_and_ders[count_samples]);
				flatoutputs.segment((count_samples)*ny,ny) = flatoutput[0];
				//flatoutputs.segment((N+1+count_samples)*ny,ny) = flatoutput[2];
				//cout<<"Flat output["<<count_samples<<"]: "<<flatoutput.transpose()<<endl;
			}
			//Tail or final state:
			{
				vector<VectorXd> flatoutput;
				(this->sys).StateAndControlsToFlatAndDerivatives(flatoutput, xs[N], us_and_ders[N-1]);//Copying us[N-1] for us[N] since us[N] does not exist
				flatoutputs.segment(N*ny,ny) = flatoutput[0];
				//flatoutputs.segment((2*N+1)*ny,ny) = flatoutput[2];
				//cout<<"Flat output["<<N<<"]: "<<flatoutput.transpose()<<endl;
			}

			//getchar();//#DEBUG
			//cout<<"Basis: "<<endl<<basis<<endl;
			//Set the output knot vector(s) as basis.pseudoinverse()*flatoutputsevaluatedatallsamples


      // Find fixed condition knot vec used to constrain solution
      VectorXd fixed_conditions;
      if(fixfinal)
      {
        fixed_conditions.resize(2*ny*(1+fixed_derivatives));
        fixed_condition_basis.resize(2*ny*(fixed_derivatives+1), (numberofknots)*ny);
      }
      else
      {
        fixed_conditions.resize(ny*(1+fixed_derivatives));
        fixed_condition_basis.resize(ny*(fixed_derivatives+1), (numberofknots)*ny);
      }

      fixed_condition_basis.block(0,0,ny*(fixed_derivatives+1), (numberofknots)*ny) = 
        CreateAllDerivativeBasis(0, tf);

      vector<VectorXd> init_conditions_vec;
      (this->sys).StateAndControlsToFlatAndDerivatives(init_conditions_vec, xs[0], us_and_ders[0]);
      assert(init_conditions_vec.size() >= 1+numberofderivatives);
      for(int i = 0; i < 1+fixed_derivatives; i++)
      {
        fixed_conditions.segment(i*ny, ny) = init_conditions_vec[i];
      }
      
      if(fixfinal)
      {
        fixed_condition_basis.block(ny*(fixed_derivatives+1),0,ny*(fixed_derivatives+1), 
          (numberofknots)*ny) = 
            CreateAllDerivativeBasis(1, tf);

        vector<VectorXd> final_conditions_vec;
        (this->sys).StateAndControlsToFlatAndDerivatives(final_conditions_vec, xs[N], us_and_ders[N-1]);
        for(int i = 0; i < 1+fixed_derivatives; i++)
        {
          fixed_conditions.segment((i+1+fixed_derivatives)*ny, ny) = final_conditions_vec[i];
        }
      }
      fixed_condition_knot_vec = 
        (fixed_condition_basis.transpose()*fixed_condition_basis).ldlt().solve(fixed_condition_basis.transpose()*fixed_conditions); 

      fixed_condition_nullspace = CalcNullspace(fixed_condition_basis);

      W.resize(2*ny*(N+1), 2*ny*(N+1));
      W.setIdentity();
      //W.block(0,0,ny*(N+1),ny*(N+1)).setZero();
      W.block(ny*(N+1),ny*(N+1),ny*(N+1),ny*(N+1)) *= 0;//.001;
      
      //HACK FOR HROTOR MPC: REMOVE THIS IN REAL LIFE!
      //for(int i = 0; i < N+1; i++)
      //{
      //  W(ny*(N+1+i)+3,ny*(N+1+i)+3) = 0;
      //}

      Eigen::MatrixXd BA = basis*fixed_condition_nullspace;
      s = (BA.transpose()*W*BA).ldlt().solve(BA.transpose()*W*(flatoutputs-basis*fixed_condition_knot_vec));
      std::cout << "s = " << s.transpose() << std::endl;
      std::cout << "size nullspace = " << fixed_condition_nullspace.rows() << "x" 
        << fixed_condition_nullspace.cols() << std::endl;
      //std::cout << "nullspace = " << std::endl << fixed_condition_nullspace << std::endl;
      std::cout << "fixed condition vec = " << fixed_conditions.transpose() << std::endl;
      std::cout << "fixed condition knot vec = " << fixed_condition_knot_vec.transpose() << std::endl;
      //std::cout << "fixed condition basis = " << std::endl << fixed_condition_basis << std::endl;

			cout<<"Error: "
        <<(basis*(fixed_condition_nullspace*s + fixed_condition_knot_vec)-flatoutputs).squaredNorm()/flatoutputs.size()<<endl;

      //cout<<"flatoutputs.head(ny) " << flatoutputs.head(ny).transpose() << endl; 
      //cout<<"flatoutputs.tail(ny) " << flatoutputs.tail(ny).transpose() << endl; 
      //(this->sys).StateAndControlsToFlat(knotsforallderivatives[0][0], xs[0], us[0]);

			//cout<<"s: "<<s.transpose()<<endl;
			//cout<<"s_temp: "<<s_temp.transpose()<<endl;
			//s.setZero();
      return true;
		}
  
  template <typename T, int nx, int nu, int np, int _ntp> 
    bool FlatOutputTparam<T, nx, nu, np, _ntp>::From(vector<double> &ts, 
                                                  vector<T> &xs, 
                                                  vector<Vectorcd> &us,
                                                  const Vectorntpd &s,
																									Vectormd *p) {
      
			assert(this->ntp == (numberofknots-(fixfinal?2:1)*(1+fixed_derivatives))*ny);
			int N = us.size();
			double tf = ts.back();//Replace this with deltat version #TODO
			vector<VectorXd> flatoutputsandderivatives(numberofderivatives+1);
			flatoutputsandderivatives_all.clear();
      flatoutputsandderivatives_all.resize(N+1);
      
      std::clock_t start = std::clock();

      // Constrain knot vector to satisfy initial (and final if specified) conditions
      Eigen::VectorXd nullspace_knot_vec = fixed_condition_nullspace*s;
      Eigen::VectorXd knot_vec = fixed_condition_knot_vec + nullspace_knot_vec;
			Createknotsforallderivatives(knot_vec);//Create all the knots before starting evaluation

			//Createknotsforallderivatives(s);//Create all the knots before starting evaluation
			//Evaluate Flat outputs from the knot inputs (s) at all the input times ts
			for(int count_ts = 0; count_ts <= N; count_ts++)
			{
				double factor_derivative = 1;
				for(int count_derivatives =0; count_derivatives <= numberofderivatives ; count_derivatives++)
				{
					if(numberofknots-count_derivatives>0)//If there are any derivatives to be evaluated then only go further
					{
						flatoutputsandderivatives[count_derivatives] = 
              (factor_derivative)*(BinomialInterp(knotsforallderivatives[count_derivatives], 
              (ts[count_ts]/tf), 0, numberofknots-count_derivatives-1));
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
				//if(count_ts == 0)
				//{
				//	T xtemp;
				//	(this->sys).FlatToStateAndControls( xtemp, us[count_ts], flatoutputsandderivatives);
				//}
        flatoutputsandderivatives_all[count_ts] = flatoutputsandderivatives;
				if(count_ts == N)
				{
					std::vector<Vectorcd> utemp;//For last point we do not care abt control
					(this->sys).FlatToStateAndControls( xs[count_ts], utemp, flatoutputsandderivatives);
				}
				else
				{
          std::vector<Vectorcd> utemp;
					(this->sys).FlatToStateAndControls( xs[count_ts], utemp, flatoutputsandderivatives);
          us[count_ts] = utemp[0];
				}
				//cout<<"us["<<count_ts<<"]: "<<us[count_ts]<<endl;
			}

			//cout<<"s: "<<s.transpose()<<endl;
                        return true;
		}
		
}

#endif
