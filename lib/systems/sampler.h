// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_SAMPLER_H
#define GCOP_SAMPLER_H

namespace gcop {
    
  /**
   * Generic sampler interface, i.e. provides a function that
   * will generate (sample) an object of type T
   *
   * Author: Marin Kobilarov
   */
  template <typename T> class Sampler {
  public:
    /**
     * Generate an object T
     * @param o object 
     * @return true if successful
     */
    virtual bool operator()(T& o) = 0;
  };
}

#endif
