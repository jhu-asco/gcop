set(sources
)

set(headers 
  shell.h
  cylinder.h
  constraint.h
)

IF (PQP_FOUND)
set(headers 
  ${headers}
  pqpdem.h
)
ENDIF (PQP_FOUND)
