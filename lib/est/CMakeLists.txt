set(headers 
  predictor.h
  corrector.h
  kalmanpredictor.h
  kalmancorrector.h
  unscentedbase.h
  unscentedpredictor.h
  unscentedcorrector.h
  gp.h
  gmm.h
  ce.h
)

set(sources 
  gp.cc
)

if(OPENCV_FOUND)
  set(headers ${headers} camera.h)
  set(sources ${sources} camera.cc)
endif(OPENCV_FOUND)

if(CERES_FOUND)
  set(headers ${headers} qrotorsystemid.h)
  set(sources ${sources} qrotorsystemid.cpp)
endif(CERES_FOUND)
