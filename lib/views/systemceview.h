// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_SYSTEMCEVIEW_H
#define GCOP_SYSTEMCEVIEW_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include "systemce.h"
#include "view.h"
#include <cmath>
#include <limits>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
 
  /**
   * SystemCe Veiw
   *
   * Authors: Marin Kobilarov
   */
  template <typename T, int n = Dynamic, int c = Dynamic, int np = Dynamic, int ntp = Dynamic> 
    class SystemCeView : public View {
    
    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, np, 1> Vectormd;

    typedef Matrix<double, n, n> Matrixnd;
    typedef Matrix<double, n, c> Matrixncd;
    typedef Matrix<double, c, n> Matrixcnd;
    typedef Matrix<double, c, c> Matrixcd;  

    typedef Matrix<double, ntp, 1> Vectortpd;  
    typedef Matrix<double, ntp, ntp> Matrixtpd;  
    
  public:
    SystemCeView(SystemCe<T, n, c, np, ntp> &ce,
                 SystemView<T, Vectorcd> &view);
    
    virtual void Render();
    
    virtual bool RenderFrame(int i);   
    
    SystemCe<T, n, c, np, ntp> &ce;
    
    SystemView<T, Vectorcd> &view;

  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int n, int c, int np, int ntp> 
    SystemCeView<T, n, c, np, ntp>::SystemCeView(SystemCe<T, n, c, np, ntp> &ce, 
                                                 SystemView<T, Vectorcd> &view) : 
    ce(ce), view(view)
  {
  }
  
  template <typename T, int n, int c, int np, int ntp> 
    void SystemCeView<T, n, c, np, ntp>::Render() {
    
    for (int j = 0; j < ce.xsas.size() && j < ce.usas.size(); ++j) {
      view.xs = &ce.xsas[j];
      view.us = &ce.usas[j]; 
      view.Render();
    }
  }

  template <typename T, int n, int c, int np, int ntp> 
    bool SystemCeView<T, n, c, np, ntp>::RenderFrame(int i) {    
    bool res = false;
    for (int j = 0; j < ce.xsas.size() && j < ce.usas.size(); ++j) {
      view.xs = &ce.xsas[j];
      view.us = &ce.usas[j]; 
      res = (view.RenderFrame(i) || res);
    }
    return res;  // return false only after all frames are rendered
  }
}

#endif
