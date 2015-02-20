// -*- coding: utf-8
// vim: set fileencoding=utf-8

// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Thomas Capricelli <orzel@freehackers.org>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Modified to find numerical difference using samples

// Author2: Gowtham

#ifndef EIGEN_SAMPLE_NUMERICAL_DIFF_H
#define EIGEN_SAMPLE_NUMERICAL_DIFF_H

#include <random>

namespace Eigen { 

/**
  * This class allows you to add a method df() to your functor, which will 
  * use numerical differentiation to compute an approximate of the
  * derivative for the functor. Of course, if you have an analytical form
  *
  * More information on
  * http://en.wikipedia.org/wiki/Numerical_differentiation
  *
  */
template<typename _Functor>
class SampleNumericalDiff : public _Functor
{
public:
    typedef _Functor Functor;
    typedef typename Functor::Scalar Scalar;
    typedef typename Functor::InputType InputType;
    typedef typename Functor::ValueType ValueType;
    typedef typename Functor::JacobianType JacobianType;

    SampleNumericalDiff(Scalar _epsfcn=0.) : Functor(), epsfcn(_epsfcn), normal_dist(0,1) {
      randgenerator.seed(370212);
      }
    SampleNumericalDiff(const Functor& f, Scalar _epsfcn=0.) : Functor(f), epsfcn(_epsfcn), normal_dist(0,1) {
      randgenerator.seed(370212);
      }

    // forward constructors
    template<typename T0>
        SampleNumericalDiff(const T0& a0) : Functor(a0), epsfcn(0) {}
    template<typename T0, typename T1>
        SampleNumericalDiff(const T0& a0, const T1& a1) : Functor(a0, a1), epsfcn(0) {}
    template<typename T0, typename T1, typename T2>
        SampleNumericalDiff(const T0& a0, const T1& a1, const T2& a2) : Functor(a0, a1, a2), epsfcn(0) {}

    enum {
        InputsAtCompileTime = Functor::InputsAtCompileTime,
        ValuesAtCompileTime = Functor::ValuesAtCompileTime
    };

    /**
      * return the number of evaluation of functor
     */
    int df(const InputType& _x, JacobianType &jac)
    {
      using std::sqrt;
      using std::abs;
      using std::cout;
      using std::endl;
      /* Local variables */
      Scalar h;
      const typename InputType::Index n = _x.size();
      int nfev= round(n*1.5);
      const Scalar eps = sqrt(((std::max)(epsfcn,NumTraits<Scalar>::epsilon() )));
      ValueType val1, val2;
      InputType x = _x;
      // TODO : we should do this only if the size is not already known
      val1.resize(Functor::values());
      val2.resize(Functor::values());
      MatrixXd dymatrix(Functor::values(), nfev);
      MatrixXd dumatrix(n, nfev);

      // initialization
      Functor::operator()(x, val1);//Mean input and value

      // Function Body
      for(int ns = 0; ns < nfev; ++ns)
      {
        for (int j = 0; j < n; ++j) {

          /*h = eps * abs(x[j]);
          if (h == 0.) {
            h = eps;
          }
          h = eps;
          if(bernoulli_dist(randgenerator))
          {
            dumatrix(j, ns) = h;
          }
          else
          {
            dumatrix(j, ns) = -h;
          }
          */
          dumatrix(j, ns) = eps*normal_dist(randgenerator);
        }//Do perturbations to the whole vector
        x = _x + dumatrix.col(ns);
        Functor::operator()(x, val2);//Evaluate at new perturbed vector
        dymatrix.col(ns) = val2 - val1;//dY
        //cout<<"dumatrix["<<ns<<"]: "<<dumatrix.col(ns).transpose()<<endl;
        //getchar();
        //cout<<"dymatrix["<<ns<<"]: "<<dymatrix.col(ns).transpose()<<endl;
        //getchar();
      }
      jac = (dumatrix*dumatrix.transpose()).ldlt().solve(dumatrix*dymatrix.transpose()).transpose();
      //cout<<"jac: "<<endl<<jac<<endl;
      cout<<endl<<"Error_predicted: "<<sqrt((dymatrix - jac*dumatrix).squaredNorm())<<endl;
      getchar();
      return nfev;
    }
private:
    Scalar epsfcn;

    SampleNumericalDiff& operator=(const SampleNumericalDiff&);

		std::default_random_engine randgenerator; ///Default random engine

		//default constructor gives p = 0.5 benoulli
		//std::bernoulli_distribution bernoulli_dist;     ///< Bernoulli generator for getting perturbations
    std::normal_distribution<double> normal_dist;///< Creates a normal distribution
};

} // end namespace Eigen

//vim: ai ts=4 sts=4 et sw=4
#endif // EIGEN_NUMERICAL_DIFF_H

