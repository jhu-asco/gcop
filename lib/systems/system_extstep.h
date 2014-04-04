#ifndef GCOP_SYSTEM_EXTSTEP_H
#define GCOP_SYSTEM_EXTSTEP_H
#include "system.h"
#include <iostream>

namespace gcop {
	template <typename T, typename Tu, int _n = Dynamic, int _c = Dynamic> 
		class System_extstep : public System<T,Tu,_n,_c>
	{
		public:

			typedef Matrix<double, _n, 1> Vectornd;
			typedef Matrix<double, _c, 1> Vectorcd;
			typedef Matrix<double, _n, _n> Matrixnd;
			typedef Matrix<double, _n, _c> Matrixncd;
			typedef Matrix<double, _c, _n> Matrixcnd;
			typedef Matrix<double, _c, _c> Matrixcd;

			typedef Matrix<double, Dynamic, 1> Vectormd;
			typedef Matrix<double, Dynamic, Dynamic> Matrixmd;
			typedef Matrix<double, _n, Dynamic> Matrixnmd;
			typedef Matrix<double, Dynamic, _n> Matrixmnd;

		private: 
			typedef boost::function<void(T &, const Tu &, double)> Func_type;
			typedef boost::function<void(void)> resetFuncType;
			Func_type  extstep;
			resetFuncType  extreset;
			//void (*extstep)(T &, const Tu &, double);// Outstate, Controls, Timestep, 

		public: 
			System_extstep(Manifold<T, _n> &X, Manifold<Tu, _c> &U, Func_type steparg, resetFuncType resetarg = NULL): System<T, Tu, _n, _c>(X,U) 
																																																					,extstep(steparg), extreset(resetarg){}

			double Step(T &xb, double t, const T &xa,
					const Tu &u, double h, const Vectormd *p = 0, 
					Matrixnd *A = 0, Matrixncd *B = 0, Matrixnmd *C = 0)
			{
				if(!extstep)
				{
					std::cout<<"external step not defined"<<std::endl;
					return 0;
				}
				extstep(xb, u, -1);// We dont need to set the timestep everytime
				return 0;
			}

			void reset()
			{
				if(!extreset)
				{
					std::cout<<"external reset not defined"<<std::endl;
					return;
				}
				extreset();
			}
	};
}

/*
	 template <typename T, typename Tu, int _n, int _c>
	 System_extstep<T,Tu, _n, _c>::System_extstep(Manifold<T, _n> &X, Manifold<Tu, _c> &U, (void)(*steparg)(const Tu &, double)):
	 System(Manifold<T, _n> &X, Manifold<Tu, _c> &U), extstep(steparg)
	 {
	 }
	 template <typename T, typename Tu, int _n, int _c> 
	 double System_extstep<T, Tu, _n, _c>::Step(T& xb, double t, const T& xa,
	 const Tu &u, double h,
	 const VectorXd *p,
	 Matrix<double, _n, _n> *A, Matrix<double, _n, _c> *B, 
	 Matrix<double, _n, Dynamic> *C) 
	 {
	 if(!extstep)
	 {
	 std::cout<<"external step not defined"<<std::endl;
	 return 0;
	 }
	 (*extstep)(xb, u, h);
	 return 0;
	 }
 */
#endif
