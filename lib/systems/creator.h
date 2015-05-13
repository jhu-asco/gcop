// This file is part of libgcop, a library for Geometric Control, Optimization, and Planning (GCOP)
//
// Copyright (C) 2004-2014 Marin Kobilarov <marin(at)jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GCOP_CREATOR_H
#define GCOP_CREATOR_H

namespace gcop {
    
  /**
   * Generic creator interface, i.e. provides a function that
   * will generate an object of type T
   *
   * Author: Marin Kobilarov
   */
  template <typename T> class Creator {
  public:
    /**
     * Generate an object T
     * @return created object
     */
    virtual T& operator()() = 0;
  };
}

#endif
