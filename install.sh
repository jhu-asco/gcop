#!/bin/bash
git submodule init #Initialize submodules
git submodule update #Clones bullet3 specific commit needed for this
cd bullet3/build3 #Go into bullet3 for building
cmake -DBUILD_SHARED_LIBS=ON -DUSE_DOUBLE_PRECISION=ON .. & make -j `nproc`
#cd ../.. & mkdir build # Go back to gcop and create a build directory
#cmake -DUSE_BULLET=ON .. & make -j `nproc`
