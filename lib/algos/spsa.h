#ifndef GCOP_SPSA_H
#define GCOP_SPSA_H

#include <Eigen/Dense>
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
 
  template <typename T, int n = Dynamic, int c = Dynamic, int _np = Dynamic> class SPSA {
    
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
		 * Implementing a Simultaneous Perturbation and Stochastic Approximation algorithm
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
    SPSA(System<T, n, c, _np> &sys, Cost<T, n, c, _np> &cost, 
         vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
         Vectormd *p = 0,          
         bool update = true);
    
    virtual ~SPSA();
    

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

    std::vector<Vectorcd> dus;      ///< control variation (N) vector

    std::vector<T> xss;             ///< states (N+1) vector
    
    std::vector<Vectorcd> uss;      ///< controls (N) vector

		std::default_random_engine randgenerator; ///Default random engine

		//default constructor gives p = 0.5 benoulli
		std::bernoulli_distribution bernoulli_dist;     ///< Bernoulli generator for getting perturbations

    int N;        ///< number of discrete trajectory segments

		int Nit;			///< Number of steps of SPSA to complete in one Iterate function

		struct Stepcoeffs{
		double a; ///< step sizes for SPSA is given ak = a/(k+1+A)^alpha
		
		double A;///< step sizes for SPSA is given ak = a/(k+1+A)^alpha

		double alpha;///< step sizes for SPSA is given ak = a/(k+1+A)^alpha
		
		double c1; ///< gradient step size is given by ck = c/(k+1)^gamma 

		double gamma; ///< gradient step size is given by ck = c/(k+1)^gamma 
		} stepc ;

    double J;    ///< optimal cost
    
    bool debug;  ///< whether to display debugging info    

		int prevcount; ///< for now use this to propagate the count forward for every iteration

  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int n, int c, int _np> 
    SPSA<T, n, c, _np>::SPSA(System<T, n, c, _np> &sys, 
                             Cost<T, n, c, _np> &cost, 
                             vector<double> &ts, 
                             vector<T> &xs, 
                             vector<Matrix<double, c, 1> > &us,
                             Matrix<double, _np, 1> *p,
                             bool update) : 
    sys(sys), cost(cost), ts(ts), xs(xs), us(us), N(us.size()), xss(xs), uss(us)
    ,Nit(200), debug(true), prevcount(0)//Choosing a, A can be done adaptively TODO
		{
			assert(N > 0);
			assert(ts.size() == N+1);
			assert(xs.size() == N+1);
			assert(us.size() == N);
			dus.resize(N);//Resizing the control variations to be same size as us

			stepc.a = 0.01;
			stepc.A = 0.1*Nit;
			stepc.c1 = 0.001;
			stepc.alpha = 0.602;
			stepc.gamma = 0.101;//Initialize the step coeffs


			if (update) {
				J = Update(xs, us);//Initial cost of the trajectory
			} else {
				J = numeric_limits<double>::max();
			}
		}
  
  template <typename T, int n, int c, int _np> 
    SPSA<T, n, c, _np>::~SPSA()
    {
    }

    
  template <typename T, int n, int c, int _np> 
    double SPSA<T, n, c, _np>::Update(vector<T> &xs, const vector<Vectorcd> &us, bool evalCost) {    
    double J = 0;
    sys.Reset(xs[0],ts[0]);//gives a chance for physics engines to reset themselves. Added Gowtham 8/2/14
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
    void SPSA<T, n, c, _np>::Iterate() {
			//Reset seed of the random generator
			randgenerator.seed(370212);
			int ussize = us[0].rows();//Number of rows in each us
			if(debug)
				cout<<"Number of rows: "<<ussize<<endl;
			float ak, ck;//step sizes
			double J1, J2;

			for(int k = 0;k < Nit;k++)
			{
				ak = stepc.a/pow((prevcount+k+1+stepc.A),(stepc.alpha));
				ck = stepc.c1/pow((prevcount+k+1),(stepc.gamma));
				if(debug)
					cout<<"Step ak, ck: "<<ak<<"\t"<<ck<<"\t"<<stepc.c1<<"\t"<<pow((prevcount+k+1),(stepc.gamma))<<endl;
				for(int count1 = 0; count1 < N; count1++)
				{
					for(int count2 = 0; count2< ussize; count2++)
					{
						if(bernoulli_dist(randgenerator))
							dus[count1][count2] = 1;
						else
							dus[count1][count2] = -1;
					}
				}//Perturbation Vector - dus
				for(int count1 = 0;count1 < N;count1++)
				{
					uss[count1] = us[count1] + ck*dus[count1];
					if(debug)
						cout<<"Uss1 #"<<count1<<": "<<uss[count1].transpose()<<endl;
				}
				if(debug)
					cout<<"Xss[0]: "<<xss[0].transpose()<<endl;
				J1 = Update(xss, uss);//Find the cost of the perturbed trajectory 1

				for(int count1 = 0;count1 < N;count1++)
				{
					uss[count1] = us[count1] - ck*dus[count1];
					if(debug)
						cout<<"Uss2 #"<<count1<<": "<<uss[count1].transpose()<<endl;
				}
				if(debug)
					cout<<"Xss[0]: "<<xss[0].transpose()<<endl;
				J2 = Update(xss, uss);//Find the cost of the perturbed trajectory 2
				//Update the control vector based on the simultaneous gradient:
				for(int count1 =0; count1 < N;count1++)
				{
					us[count1] = us[count1] - ((ak/ck)*(J1 - J2)/2)*dus[count1].cwiseInverse();//This is ok only for the special case of bernoulli distribution with dus being +1 or -1 since element wise inverse gives the same values back again
				}
			}
			prevcount += Nit;
			//Terminal xs and J:
			J = Update(xs,us);
		}
}

#endif
