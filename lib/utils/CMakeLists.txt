set(sources
  utils.cc
  se2.cc
  so3.cc
  se3.cc
  params.cc
  urdflink.cpp
  urdfjoint.cpp
#  normal.cc
  )

set(headers
  utils.h
  group.h
  se2.h
  so3.h
  se3.h
	urdflink.h
	urdfjoint.h
	pose.h
	color.h
  urdfmodel.h
  function.h
  params.h
  normal.h
  bulletworld.h
  samplenumericaldiff.h
  load_eigen_matrix.h
  )

if (casadi_FOUND)
  set(headers ${headers}
        casadi_fully_connected_layer.h
        eigen_casadi_conversions.h
        fully_connected_layer.h)
  set(sources ${sources}
      fully_connected_layer.cpp)
endif (casadi_FOUND)
