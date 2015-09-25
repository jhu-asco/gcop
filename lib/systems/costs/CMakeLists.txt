set(sources
#    body2dcost.cc
    gunicyclecost.cc
    body2dtrackcost.cc
    body3dtrackcost.cc
    bacost.cc
    body2dslamcost.cc
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
    multilscost.h
    constraintcost.h
    body2dtrackcost.h
    body3dtrackcost.h
    kinbody3dtrackcost.h
    kinrccarpathcost.h
    kinbodyprojtrackcost.h
    bacost.h
    body2dslamcost.h
#airplanecost.h
)

IF (PQP_FOUND)
set(headers
  ${headers}
  pqpdemcost.h
)
ENDIF (PQP_FOUND)
