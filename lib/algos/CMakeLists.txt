set(headers
      ddp.h
      sddp.h
      docp.h
      doep.h
      gndocp.h
      gndoep.h
      deterministicgndoep.h
      pddp.h
      ba.h 
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
      body2dslam.cc
      body2dgraph.cc
      posegraph2d.cc
      mbscontroller.cc
      point3dcontroller.cc
      gavoidcontroller.cc
)

IF (CERES_FOUND)
set(headers 
  ${headers}
  dynvisins.h
)
set(sources 
  ${sources}
  dynvisins.cc
)
ENDIF (CERES_FOUND)


IF (DSL_FOUND AND PQP_FOUND)
set(headers 
  ${headers}
  body3ddemcontroller.h
)
ENDIF (DSL_FOUND AND PQP_FOUND)

IF (CERES_FOUND)
set(headers
    ${headers}
    dynvisins.h
    )
set (sources
    ${sources}
    dynvisins.cc
    )
ENDIF (CERES_FOUND)
