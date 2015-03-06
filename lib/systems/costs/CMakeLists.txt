set(sources
    body2dcost.cc
    gunicyclecost.cc
    body2dtrackcost.cc
    body3dtrackcost.cc
#airplanecost.cc
   )

set(headers
    body2dcost.h
    body3dcost.h
    kinbody3dcost.h
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
    body2dtrackcost.h
    body3dtrackcost.h
    kinbody3dtrackcost.h
    kinbodyprojtrackcost.h
#airplanecost.h
)

IF (PQP_FOUND)
set(headers
  ${headers}
  pqpdemcost.h
)
ENDIF (PQP_FOUND)
