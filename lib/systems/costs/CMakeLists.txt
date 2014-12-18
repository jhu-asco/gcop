set(sources
    body2dcost.cc
    gunicyclecost.cc
#airplanecost.cc
   )

set(headers
    body2dcost.h
    body3dcost.h
    cost.h
    gunicyclecost.h
    lqcost.h
    lqsensorcost.h
    lscost.h
    lssensorcost.h
    mbslqcost.h
    rnlqcost.h
    rnlqsensorcost.h
    sensorcost.h
    multicost.h
    constraintcost.h
#airplanecost.h
)

IF (PQP_FOUND)
set(headers
  ${headers}
  pqpdemcost.h
)
ENDIF (PQP_FOUND)
