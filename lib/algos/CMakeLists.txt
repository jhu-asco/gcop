set(headers
      ddp.h
      sddp.h
      docp.h
      doep.h
      gndocp.h
      gndoep.h
      gndoepv2.h
      pddp.h
      ba.h
      body2dslam.h
      bacost.h
      body2dslamcost.h
      body2dgraph.h
      posegraph2d.h
      controller.h
      mbscontroller.h
      point3dcontroller.h
      body3dcontroller.h
      gavoidcontroller.h
      body3davoidcontroller.h
      systemce.h
      aspsa.h
      spsa.h
)

set(sources
      ba.cc
      bacost.cc
      body2dslam.cc
      body2dslamcost.cc
      body2dgraph.cc
      posegraph2d.cc
      mbscontroller.cc
      point3dcontroller.cc
      gavoidcontroller.cc
)

IF (DSL_FOUND AND PQP_FOUND)
set(headers 
  ${headers}
  body3ddemcontroller.h
)
ENDIF (DSL_FOUND AND PQP_FOUND)