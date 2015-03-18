#ifndef GCOP_SYSTEM_EXTSTEP_H
#define GCOP_SYSTEM_EXTSTEP_H
#include "system.h"
#include <iostream>

namespace gcop {
  template <typename T, 
    int _nx = Dynamic, 
    int _nu = Dynamic,
    int _np = Dynamic> 
    class System_extstep : public System<T,_nx,_nu,_np>
    {
      public:
        typedef Matrix<double, _nx, 1> Vectornd;
        typedef Matrix<double, _nu, 1> Vectorcd;
        typedef Matrix<double, _np, 1> Vectormd;

        typedef Matrix<double, _nx, _nx> Matrixnd;
        typedef Matrix<double, _nx, _nu> Matrixncd;
        typedef Matrix<double, _nu, _nx> Matrixcnd;
        typedef Matrix<double, _nu, _nu> Matrixcd;

        typedef Matrix<double, _np, _np> Matrixmd;
        typedef Matrix<double, _nx, _np> Matrixnmd;
        typedef Matrix<double, _np, _nx> Matrixmnd;

      private: 
        typedef boost::function<void(T &, const Vectorcd &, double)> Func_type;
        typedef boost::function<void(const T &)> resetFuncType;
        Func_type  extstep;
        resetFuncType  extreset;
        //void (*extstep)(T &, const Tu &, double);// Outstate, Controls, Timestep, 

      public: 
        //np is not added but can be added if needed
        System_extstep(Manifold<T, _nx> &X, int nu, Func_type steparg, resetFuncType resetarg = NULL): System<T,_nx,_nu,_np>(X,nu) 
                                                                                                       ,extstep(steparg), extreset(resetarg){}

        double Step(T &xb, double t, const T &xa,
            const Vectorcd &u, double h, const Vectormd *p = 0, 
            Matrixnd *A = 0, Matrixncd *B = 0, Matrixnmd *C = 0, bool reset=false)
        {
          if(!extstep || !extreset)
          {
            std::cout<<"external step or reset not defined"<<std::endl;
            return 0;
          }
          extreset(xa);//Reset to initial state
          extstep(xb, u, h);
          //We do not update internal state in finite difference functions
          return 0;
        }

        double Step_internaloutput(const Vectorcd &u, double h,
            const Vectormd *p,
            Matrixnd *A, Matrix<double, _nx, _nu> *B, 
            Matrix<double, _nx, _np> *C) 
        {
          if(!extstep)
          {
            std::cout<<"External step not defined"<<std::endl;
          }
          extstep(this->x, u, h);
          this->t = (this->t + h);
          return 0;
        }

        double Step_internalinput(T& xb,
            const Vectorcd &u, double h,
            const Vectormd *p,
            Matrixnd *A, Matrix<double, _nx, _nu> *B, 
            Matrix<double, _nx, _np> *C) 
        {
          extstep(xb, u, h);
          this->x = xb;
          this->t = (this->t + h);
          return 0;
        }

        bool reset(const T& x, double t = 0)
        {
          extreset(x);
          this->x = x;
          this->t = t;
        }
    };
}
#endif
