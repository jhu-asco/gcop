set(headers 
  predictor.h
  corrector.h
  kalmanpredictor.h
  kalmancorrector.h
  unscentedbase.h
  unscentedpredictor.h
  unscentedcorrector.h
  quat.h
  gp.h
  gmm.h
  ce.h
  imu.h
  imumanifold.h
  sensor.h
  imusensor.h
  ins.h
  insmanifold.h
  insimu.h
  insgps.h
  point3d.h
  point3dmanifold.h
  point3dgps.h
)

if(OPENCV_FOUND)
  set(headers ${headers} camera.h)
endif(OPENCV_FOUND)


set(sources 
  quat.cc
  gp.cc
  imu.cc
  imumanifold.cc
  ins.cc
  insmanifold.cc
  point3d.cc
  point3dmanifold.cc
)

if(OPENCV_FOUND)
  set(sources ${sources} camera.cc)
endif(OPENCV_FOUND)
