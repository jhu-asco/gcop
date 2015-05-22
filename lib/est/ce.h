// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_CE_H
#define GCOP_CE_H

#include <vector>
#include "gmm.h"
#include <limits>

namespace gcop {
  
  using namespace Eigen;
  using namespace std;
  
  /**
   * Cross-entropy optimization method. One iteration of the algorithm typically includes the following steps:
   *
   * 0. Initialize the GMM distribution, select number of samples N, and quantile rho
   *
   * 1. for 1:N, z = Sample(), c = Cost(z), AddSample(z,c), end
   *
   * 2. Select() selects the top rho-quantile by sorting the samples using their costs
   *
   * 3. Fit() estimates the distribution
   *
   * 4. Best() gives the best sample or just take the first sample in the list of samples as the solution
   *
   * 5. Reset() and goto 1.
   *
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  template <int _n = Dynamic> 
    class Ce {
  public:  
  typedef Matrix<double, _n, 1> Vectornd;
  typedef Matrix<double, _n, _n> Matrixnd;
  
  /**
   * Initialize n-dimensional CE with k modes
   * @param n parameter space dimension
   * @param k number of mixture components
   * @param S n-n matrix of random noise to be added to the covariance after GMM estimation; acts both as regularizing factor and also for improving local exploration. This should typically be set to a small value, e.g. a diagonal matrix with entries = 0.01
   */    
  Ce(int n = _n, int k = 1, const Matrixnd *S = 0);
  
  virtual ~Ce();
  
  /**
   * Clear samples
   */
  void Reset();
  
  /**
   * Add a sample to list of samples
   * @param z sample
   * @param c cost
   */
  virtual void AddSample(const Vectornd &z, double c);
  
  /**
   * Select the top quantile. Resizes samples to size N*rho containing the best samples
   */
  void Select();
  
  /**
   * Fit GMM to data
   */
  bool Fit();
  
  /**
   * Draw a sample from the GMM
   * @param z vector to be sampled
   * @return likelihood
   */
  virtual double Sample(Vectornd &z);
  
  /**
   * Get the best sample. Equivalent to accessing zps[0].first. This assumes
   * Select() was called to sort the samples.
   * @return the sample with lowest cost
   */
  const Vectornd& Best();
  
  int n;          ///< parameter space dimension
      
  Gmm<_n> gmm;    ///< underlying Gaussian-mixture-model

  Matrixnd S;     ///< extra noise

  vector<pair<Vectornd, double> > zps;  ///< N samples (pair of vector and its likelihood)

  vector<double> cs;    ///< costs of the N samples

  vector<pair<Vectornd, double> > zpes;  ///< elite samples (pair of vector and its likelihood)

  double rho;     ///< quantile: should be b/n .01 and .1 (default is .1)

  double alpha;   ///< update factor: parameter is updated according to v=a*v_new + (1-a)*v_old , i.e. alpha is used to "smooth" the update using the old value (default is .9)

  bool mras;      ///< whether to use MRAS version of algorithm, i.e. instead of using rare-event volume estimation, tilt the density optimally using the Gibbs density ~exp(-b*J(z))

  double b;       ///< Gibbs factor (only applicable if mras = true)

  bool bAuto;    ///< let b be determined automatically (true by default)

  bool inc;      ///< incremental version of mras

  Vectornd zmin; ///< current best

  double Jmin;    ///< minimum cost

  double Jmax;    ///< current Jmax (discard any samples with costs higher than Jmax)

  double Jave;    ///< average cost over all samples

  };
  
  template <int _n>
    Ce<_n>::Ce(int n,
               int k,
               const Matrixnd *S) :
    n(n),
    gmm(n, k),
    rho(.1),
    alpha(.9),
    mras(false),
    b(1),
    bAuto(true),
    inc(false),
    Jmin(std::numeric_limits<double>::max()),
    Jmax(0),
    Jave(std::numeric_limits<double>::max())
      {
        if (_n == Dynamic) {
          this->S.resize(n, n);
          this->zmin.resize(n);
        } else {
          assert(_n == n);
        }
        
        if (S)
          this->S = *S;
        else
          this->S.setZero();
      }
  
  template <int _n>
    Ce<_n>::~Ce()
    {
    }


  template <int _n>
    void Ce<_n>::Reset()
    {
      zps.clear();
      cs.clear();
      Jmin = std::numeric_limits<double>::max();
      Jmax = 0;
    }

  template <int _n>
    void Ce<_n>::AddSample(const Vectornd &z, double c) 
    { 
      int Ns = zps.size();
      if (!Ns) {
        Jave = c;
        Jmax = c;
      } else {
        Jave = (Jave*Ns + c)/(Ns + 1);
      }
      zps.push_back(make_pair(z, c));
      cs.push_back(c);

      if (Jmin > c) {
        zmin = z;
        Jmin = c;
      }

      if (Jmax < c) {
        Jmax = c;
      }      
    }

  template <int _n>
    bool zpSort(const pair<Matrix<double, _n, 1>, double> &zpa, const pair<Matrix<double, _n, 1>, double> &zpb)
    {
      return zpa.second < zpb.second;
    }

  template <int _n>
    void Ce<_n>::Select()
    {
      if (mras) {
        // set b to minimum cost
        if (bAuto) {
          b = 10/Jmin;
          /*
          b = std::numeric_limits<double>::max();
          for (int j = 0; j < cs.size(); ++j) {xo
            if (b > cs[j])
              b = cs[j];
          }
          */
        }
        return;
      }
  
      int N = zps.size();
      int Ne = MIN((int)ceil(N*rho), N);
      zpes.clear();
  
      if (Ne < N) {
        std::sort(zps.begin(), zps.end(), zpSort<_n>);
        std::sort(cs.begin(), cs.end());   // this is redundant but is kept for consistency
        for (int i = 0; i < Ne; ++i)
          zpes.push_back(zps[i]);
        // std::copy(zpes.begin(), zps.begin(), zps.begin() + Ne);

        //        zps.resize(Ne);
        //        cs.resize(Ne);
      }
    }


  template <int _n>
    bool Ce<_n>::Fit()
    {       
      int N = zps.size();
  
      if (mras) {
        double cn = 0;
        for (int j = 0; j < N; ++j) {
          pair<Vectornd, double> &zp = zps[j];
          //          zp.second = exp(-b*(cs[j] - Jmin)/(Jmax - Jmin));     // pdf
          zp.second = exp(-b*cs[j]);
          cn += zp.second;               // normalizer
        }
        assert(cn > 0);
        for (int j = 0; j < N; ++j) 
          zps[j].second /= cn;               // normalize
    
        gmm.Fit(zps, alpha, 50, &S);
      } else {
        int Ne = zpes.size();
        for (int j = 0; j < Ne; ++j) 
          zpes[j].second = 1.0/Ne;             // probability
    
        gmm.Fit(zpes, alpha, 50, &S);
      }

      return gmm.Update();
    }

  template <int _n>
    double Ce<_n>::Sample(Vectornd &z) 
    {
      return gmm.Sample(z);
    }

  template <int _n>
    const Matrix<double, _n, 1>& Ce<_n>::Best() 
    {
      return zmin;
    }
}


#endif
