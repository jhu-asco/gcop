set(sources
)

set(headers 
  shell.h
  cylinder.h
)

IF (PQP_FOUND)
set(headers 
  ${headers}
  pqpdem.h
)
ENDIF (PQP_FOUND)