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

if(OPENCV_FOUND)
  set(headers ${headers} camera.h)
endif(OPENCV_FOUND)


set(sources 
  gp.cc
)

if(OPENCV_FOUND)
  set(sources ${sources} camera.cc)
endif(OPENCV_FOUND)
