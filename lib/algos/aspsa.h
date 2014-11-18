#ifndef GCOP_ASPSA_H
#define GCOP_ASPSA_H

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include "system.h"
#include "cost.h"
#include <cmath>
#include "rn.h"
#include <limits>
#include <random>

namespace gcop {

	using namespace std;
	using namespace Eigen;

	template <typename T, int n = Dynamic, int c = Dynamic, int _np = Dynamic> class ASPSA {

		typedef Matrix<double, n, 1> Vectornd;
		typedef Matrix<double, c, 1> Vectorcd;
		typedef Matrix<double, _np, 1> Vectormd;
		typedef Matrix<double, n, n> Matrixnd;
		typedef Matrix<double, n, c> Matrixncd;
		typedef Matrix<double, c, n> Matrixcnd;
		typedef Matrix<double, c, c> Matrixcd;  

                typedef Matrix<double, _np, _np> Matrixmd;
                typedef Matrix<double, n, _np> Matrixnmd;
                typedef Matrix<double, _np, n> Matrixmnd;


		public:
		/**
		 * Implementing Adaptive Simultaneous Perturbation and Stochastic Approximation algorithm
		 * in GCOP library. The gains needed for the algorithms should be chosen wisely using
		 * Section 7.5.2 of Introduction to Stochastic Search and Optimization (ISSO) book.
		 *
		 * This method minimizes the cost of a trajectory using the input parameters as us
		 * The times ts must be given, the initial state xs[0] must be set, and
		 * the controls us will be used as an initial guess for the optimization.
		 *
		 * After initialization, every call to Iterate() will optimize the 
		 * controls us and states xs and modify them accordingly. Problems involving
		 * time-optimization will also modify the sequence of times ts.
		 * 
		 * @param sys system
		 * @param cost cost
		 * @param ts (N+1) sequence of discrete times
		 * @param xs (N+1) sequence of discrete states
		 * @param us (N) sequence of control inputs
		 * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
		 *               This is necessary only if xs was not already generated from us.
		 */
		ASPSA(System<T, n, c, _np> &sys, Cost<T, n, c, _np> &cost, 
                      vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
                      Vectormd *p = 0, 
                      bool update = true);
                
		virtual ~ASPSA();


		/**
		 * Perform one SYSTEMCE iteration. 
		 */
		void Iterate();


		/**
		 * Generate a full trajectory (xs and us) from a parameter z and (optionally) return its cost
		 * @param xs trajectory
		 * @param us controls
		 * @param evalCost whether to compute and return the trajectory cost
		 */
		double Update(std::vector<T> &xs, const std::vector<Vectorcd> &us, bool evalCost = true);


		System<T, n, c, _np> &sys;    ///< dynamical system

		Cost<T, n, c, _np> &cost;     ///< given cost function

		std::vector<double> &ts; ///< times (N+1) vector

		std::vector<T> &xs;      ///< states (N+1) vector

		std::vector<Vectorcd> &us;      ///< controls (N) vector

		std::vector<Vectorcd> dus1;      ///< control variation (N) vector
		
		std::vector<Vectorcd> ustemp;      ///< control variation (N) vector

		std::vector<Vectorcd> dus2;      ///< control variation (N) vector

		std::vector<T> xss;             ///< states (N+1) vector

		std::vector<Vectorcd> uss1;      ///< controls (N) vector perturb1

		std::vector<Vectorcd> uss2;      ///< controls (N) vector perturb2
		
		std::vector<Vectorcd> ghat;      ///< Estimated gradient

		std::vector<Vectorcd> ghatavg;      ///< Estimated avg gradient

		std::vector<Vectorcd> hessian_currest;      ///< Current iterate hessian estimate

		std::vector<Vectorcd> hessian_currestavg;      ///< Current iterate hessian estimate after averaging

		std::vector<Vectorcd> hessian_est;      ///< Estimating Hessian

		std::vector<Vectorcd> hessian_estreg;      ///< Adding Regularizing term

		std::vector<Vectorcd> deltaG;      ///< Difference in Gradient

		std::default_random_engine randgenerator; ///Default random engine

		//default constructor gives p = 0.5 benoulli
		std::bernoulli_distribution bernoulli_dist;     ///< Bernoulli generator for getting perturbations

		int N;        ///< number of discrete trajectory segments

		int Nit;			///< Number of steps of ASPSA to complete in one Iterate function

		float toleranceJ; ///< Tolerance for how high J can go from the current value to reject some 

		float toleranceHessianmag; ///< Tolerance for how high J can go from the current value to reject some 

		struct Stepcoeffs{
			double a; ///< step sizes for SPSA is given ak = a/(k+1+A)^alpha

			double a1; ///< step sizes for ASPSA is given ak = a/(k+1+A)^alpha

			double A;///< step sizes for SPSA is given ak = a/(k+1+A)^alpha

			double A1;///< step sizes for SPSA is given ak = a/(k+1+A)^alpha

			double alpha;///< step sizes for SPSA is given ak = a/(k+1+A)^alpha

			double alpha1;///< step sizes for ASPSA is given ak = a/(k+1+A)^alpha

			double c1; ///< gradient step size is given by ck = c/(k+1)^gamma  for SPSA

			double c11; ///< gradient step size is given by ck = c/(k+1)^gamma  for 2SPSA

			double c2; ///< gradient step size is given by cktilde = ctilde/(k+1)^gamma  secondary perturbation 2SPSA

			double gamma; ///< gradient step size is given by ck = c/(k+1)^gamma  for SPSA

			double gamma1; ///< gradient step size is given by ck = c/(k+1)^gamma 

			double Navg;   ///< Averaging the Gradient and Hessian
		} stepc ;

		double J;    ///< optimal cost

		bool debug;  ///< whether to display debugging info    

		int prevcount; ///< for now use this to propagate the count forward for every iteration

	};

	using namespace std;
	using namespace Eigen;

	template <typename T, int n, int c, int _np> 
          ASPSA<T, n, c, _np>::ASPSA(System<T, n, c, _np> &sys, 
				Cost<T, n, c, _np> &cost, 
				vector<double> &ts, 
				vector<T> &xs, 
				vector<Matrix<double, c, 1> > &us,
                                Matrix<double, _np, 1> *p,                                
				bool update) : 
          sys(sys), cost(cost), ts(ts), xs(xs), us(us), dus1(us), dus2(us), N(us.size()), xss(xs), uss1(us), uss2(us), ustemp(us)
			,Nit(200), debug(true), prevcount(0), toleranceJ(0.0001), toleranceHessianmag(0.1)//Choosing a, A can be done adaptively TODO
	{
		assert(N > 0);
		assert(ts.size() == N+1);
		assert(xs.size() == N+1);
		assert(us.size() == N);
		//dus.resize(us.size());//Resizing the control variations to be same size as us

		stepc.a = 0.01;
		stepc.a1 = 1;
		stepc.A1 = 0;
		stepc.alpha1 = 1;
		stepc.gamma1 = 1.0/6.0;
		stepc.A = 0.1*Nit;
		stepc.c1 = 0.001;
		stepc.c11 = 2*0.001;
		stepc.c2 = 2*0.001;
		stepc.alpha = 0.602;
		stepc.gamma = 0.101;//Initialize the step coeffs
		stepc.Navg = 4;

		hessian_est.resize(N);//Resizing Hessian to be diagonal matrix
		hessian_currest.resize(N);
		hessian_currestavg.resize(N);
		hessian_estreg.resize(N);
		deltaG.resize(N);
		ghat.resize(N);
		ghatavg.resize(N);
		for(int count1 = 0;count1 < N;count1++)
		{
			hessian_est[count1] = 5000*Vectorcd::Ones();//Initialize hessian_est
		}

		if (update) {
			J = Update(xs, us);//Initial cost of the trajectory
		} else {
			J = numeric_limits<double>::max();
		}
	}

	template <typename T, int n, int c, int _np> 
          ASPSA<T, n, c, _np>::~ASPSA()
          {
          }


	template <typename T, int n, int c, int _np> 
          double ASPSA<T, n, c, _np>::Update(vector<T> &xs, const vector<Vectorcd> &us, bool evalCost) {    
          double J = 0;
          sys.reset(xs[0],ts[0]);//gives a chance for physics engines to reset themselves. Added Gowtham 8/2/14
          for (int k = 0; k < N; ++k) {
            double h = ts[k+1] - ts[k];
            sys.Step(xs[k+1], us[k], h);
            if (evalCost) 
              J += cost.L(ts[k], xs[k], us[k], h);
            
          }
          if (evalCost)
            J += cost.L(ts[N], xs[N], us[N-1], 0);
          return J;
        }
        
	template <typename T, int n, int c, int _np> 
          void ASPSA<T, n, c, _np>::Iterate() {
          //Reset seed of the random generator
          randgenerator.seed(370212);
          int ussize = us[0].rows();//Number of rows in each us
          if(debug)
            cout<<"Number of rows: "<<ussize<<endl;
          float ak, ck, cktilde;//step sizes for 1SPSA and 2SPSA
          float J1, J2, J3, J4, Jtemp;
          
          //1SPSA for Nit steps
          for(int k = 0;k < 40;k ++)
            {
              ak = stepc.a/pow((prevcount+k+1+stepc.A),(stepc.alpha));
              ck = stepc.c1/pow((prevcount+k+1),(stepc.gamma));
              if(debug)
                cout<<"Step ak, ck: "<<ak<<"\t"<<ck<<"\t"<<stepc.c1<<"\t"<<pow((prevcount+k+1),(stepc.gamma))<<endl;
              for(int count1 = 0;count1 < N;count1++)
                {
                  for(int count2 = 0; count2< ussize; count2++)
                    {
                      if(bernoulli_dist(randgenerator))
                        dus1[count1][count2] = 1;
                      else
                        dus1[count1][count2] = -1;
                    }//Perturbation Vector - dus
                  uss1[count1] = us[count1] + ck*dus1[count1];
                  uss2[count1] = us[count1] - ck*dus1[count1];
                  /*
                    if(debug)
                    {
                    cout<<"Uss1 #"<<count1<<": "<<uss1[count1].transpose()<<endl;
                    cout<<"Uss2 #"<<count1<<": "<<uss2[count1].transpose()<<endl;
                    }
                  */
                }//Perturbed once thetahatk + ck deltak
              J1 = Update(xss, uss1);//Find the cost of the perturbed trajectory 1
              J2 = Update(xss, uss2);//Find the cost of the perturbed trajectory 2
              
              
              
              //Update the control vector based on the simultaneous gradient:
              for(int count1 =0; count1 < N;count1++)
                {
                  ghat[count1] = ((J1 - J2)/(2*ck))*(dus1[count1].cwiseInverse());
                  us[count1] = us[count1] - ak*ghat[count1];
                }
              if(debug)
                {
                  J = Update(xs,us);
                  cout<<"J #"<<k<<"\t: "<<J<<endl;
                }
            }
          J = Update(xs,us);
          cout<<"J "<<"\t: "<<J<<endl;
          getchar();
          
          prevcount += 40;//Adding the number of iterations done till now
          
          //2SPSA Begin:
          for(int k = 0;k < Nit;k++)
            {
              ak = stepc.a1/pow((prevcount+k+1+stepc.A1),(stepc.alpha1));
              ck = stepc.c11/pow((prevcount+k+1),(stepc.gamma1));
              cktilde = stepc.c2/pow((prevcount+k+1),(stepc.gamma1));
              if(debug)
                cout<<"Step ak, ck, cktilde: "<<ak<<"\t"<<ck<<"\t"<<cktilde<<endl;
              //Averagind Gradient and Hessian information
              
              for(int count1 = 0;count1 < N;count1++)
                {
                  hessian_currestavg[count1].setZero();//Initialize hessian avg to 0
                  ghatavg[count1].setZero();//initialize ghatavg
                }
              
              for(int mavg = 1;mavg <= stepc.Navg;mavg++)
                {
                  for(int count1 = 0;count1 < N;count1++)
                    {
                      for(int count2 = 0; count2< ussize; count2++)
                        {
                          if(bernoulli_dist(randgenerator))
                            dus1[count1][count2] = 1;
                          else
                            dus1[count1][count2] = -1;
                        }//Perturbation Vector - delta
                      uss1[count1] = us[count1] + ck*dus1[count1];
                      uss2[count1] = us[count1] - ck*dus1[count1];
                      /*
                        if(debug)
                        {
                        cout<<"Uss1 #"<<count1<<": "<<uss1[count1].transpose()<<endl;
                        cout<<"Uss2 #"<<count1<<": "<<uss2[count1].transpose()<<endl;
                        }
                      */
                    }//Perturbed once thetahatk +/- ck deltak
                  J1 = Update(xss, uss1);//Find the cost of the perturbed trajectory 1
                  J2 = Update(xss, uss2);//Find the cost of the perturbed trajectory 2
                  
                  
                  //Second perturbation to find the gradient at thethatk + ck deltak
                  for(int count1 = 0;count1 < N;count1++)
                    {
                      for(int count2 = 0; count2< ussize; count2++)
                        {
                          if(bernoulli_dist(randgenerator))
                            dus2[count1][count2] = 1;
                          else
                            dus2[count1][count2] = -1;
                        }//Perturbation Vector - deltatilde
                      uss1[count1] += cktilde*dus2[count1];
                      uss2[count1] += cktilde*dus2[count1];
                      /*if(debug)
                        {
                        cout<<"Uss1 #"<<count1<<": "<<uss1[count1].transpose()<<endl;
                        cout<<"Uss2 #"<<count1<<": "<<uss2[count1].transpose()<<endl;
							}
                      */
                    }
                  J3 = Update(xss, uss1);//Find the cost of the additional perturbations
                  J4 = Update(xss, uss2);
                  //Estimate the gradients:
                  //Gk1(thetahatk + ck deltak) =  (y(thethatk + ckdeltak + cktilde deltak) - y(thetahatk + ck deltak))/cktilde * deltak inverse
                  //float minhessianval = 1e4, minhesscount1 = 0;
                  //Finding average Hessian and Average Gradient
                  for(int count1 =0; count1 < N;count1++)
                    {
                      deltaG[count1] = ((J3 - J1 - (J4 - J2))/cktilde)*(dus2[count1].cwiseInverse());
                      hessian_currest[count1] = (1/(2*ck))*(deltaG[count1].cwiseProduct(dus1[count1].cwiseInverse()));//Considering only diagonal elements
                      /*if(debug
                        {
                        cout<<"hesscount2 #"<<count1<<":\n "<<hessian_currest[count1]<<endl;
                        }
                      */
                      //hessian_currest[count1] = 0.5*(hessian_currest[count1] + hessian_currest[count1].transpose());
                      
                      //Computing the averages:
                      ghat[count1] = ((J1 - J2)/(2*ck))*(dus1[count1].cwiseInverse());
                      ghatavg[count1] = (float(mavg-1)/float(mavg))*ghatavg[count1] + (1.0/float(mavg))*ghat[count1];//Estimating avg of the gradient
                      hessian_currestavg[count1] = (float(mavg-1)/float(mavg))*hessian_currestavg[count1] + (1.0/float(mavg))*hessian_currest[count1];//Estimating hessian avg 
                      /*	if(debug)
                                {
                                //cout<<"us[#"<<count1<<"]: "<<us[count1].transpose()<<endl;
                                //cout<<"hesscount #"<<count1<<":\n "<<hessian_estreg[count1]<<endl;
                                cout<<"hesscount2 #"<<count1<<":\n "<<hessian_currest[count1]<<endl;
                                }*/
                    }
                  if(debug)
                    {
                      cout<<"J1,J2,J3,J4:"<<J1<<"\t"<<J2<<"\t"<<J3<<"\t"<<J4<<"\t"<<endl;
                      cout<<"Magof Hessian"<<(J3 - J1 - (J4 - J2))/(2*cktilde*ck)<<"\t First Gradient: "<<(J3 - J1)/(cktilde)<<"\t Second Gradient: "<<(J4 - J2)/(cktilde)<<endl;
                      cout<<"Magof Gradient"<<((J1 - J2)/(2*ck))<<endl;
                      getchar();
                    }
                  //					if(debug)
                  //getchar();//DEBUG STUFF
                }//Averaging Done
              
              float magofHessian = abs((J3 - J1 - (J4 - J2))/(2*cktilde*ck));//Just checking one of it instead of the actual average TODO
              /*if(magofHessian >toleranceHessianmag)
					{
					cerr<<"Magnitude of Hessian Big: "<<magofHessian<<endl;
					for(int count1 =0; count1 < N;count1++)
					{
					ustemp[count1] = us[count1] - ak*ghatavg[count1];//Can add constraints here TODO
					}
					}
					else
					{
              */
              //Update the control vector based on the simultaneous gradient:
              for(int count1 =0; count1 < N;count1++)
                {
                  hessian_est[count1] = (float(k+1)/float(k+2))*hessian_est[count1] + (1/float(k+2))*hessian_currestavg[count1];//Find the Current Hessian estimate using iterative process This is slightly modified to use the initial estimate 
                  //hessian_estreg[count1] = hessian_est[count1].cwiseAbs() + (0.001/(k+1)) * Vectorcd::Ones();//Regulated
                  hessian_estreg[count1] = hessian_est[count1].cwiseMax(toleranceHessianmag*Vectorcd::Ones());//Regulated
                  //hessian_estreg[count1] = Vectorcd::Ones();
                  ustemp[count1] = us[count1] - ak*(hessian_estreg[count1].cwiseInverse()).cwiseProduct(ghatavg[count1]);//Can add constraints here TODO
                  if(debug)
                    {
                      //cout<<"us[#"<<count1<<"]: "<<us[count1].transpose()<<endl;
                      //cout<<"Hess eig values: "<<hessian_estreg[count1].eigenvalues().transpose()<<endl;
                      cout<<"hesscountreg #"<<count1<<": "<<hessian_estreg[count1].transpose()<<endl;
                      cout<<"hesscountavgest #"<<count1<<": "<<hessian_currestavg[count1].transpose()<<endl;
                      //.part<Eigen::LowerTriangular>() 
                      //cout<<"hesscount2 #"<<count1<<":\n "<<hessian_currest[count1]<<endl;
                    }
                }
              //}
              Jtemp = Update(xs,ustemp);//Evaluate the new cost
              if(Jtemp < (J + toleranceJ))
                {
                  J = Jtemp;
                  us = ustemp;	
                }//Else neglect
              else
                {
                  cerr<<"Neglecting Jtemp, J: "<<Jtemp<<"\t"<<J<<endl;
                }
              if(debug)
                {
                  cout<<"J Jtemp J1,J2,J3,J4 #"<<k<<"\t: "<<J<<"\t"<<Jtemp<<"\t"<<J1<<"\t"<<J2<<"\t"<<J3<<"\t"<<J4<<endl;
                  getchar();
                }
            }
          prevcount += Nit;
          //Terminal xs and J:
          J = Update(xs,us);
        }
}

#endif
