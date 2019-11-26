set(sources
  disk.cc
  sphere.cc
)

set(headers 
  shell.h
  cylinder.h
  disk.h
  sphere.h
  diskconstraint.h
  sphereconstraint.h
  constraint.h
  groundplane.h
  yawvelfollow.h
)

IF (PQP_FOUND)
set(headers 
  ${headers}
  pqpdem.h
)
ENDIF (PQP_FOUND)
