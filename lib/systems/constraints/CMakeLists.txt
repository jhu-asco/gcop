set(sources
  disk.cc
)

set(headers 
  shell.h
  cylinder.h
  disk.h
  diskconstraint.h
  constraint.h
)

IF (PQP_FOUND)
set(headers 
  ${headers}
  pqpdem.h
)
ENDIF (PQP_FOUND)
