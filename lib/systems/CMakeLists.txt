set(sources
    particle2d.cc
    unicycle.cc
    car.cc
    rccar.cc
    gunicycle.cc
    qrotor.cc
    hrotor.cc
    heli.cc
    kinbody2d.cc
    body2d.cc
    chain.cc
    chain1.cc
    airm.cc
    airbot.cc
    airbase.cc
    mbs.cc
    snake.cc
    joint.cc
    dem.cc
    joint.cc
    body2dforce.cc
    arm.cc
    urdf_parser.cpp
    ins.cc
    point3d.cc
    imu.cc
)

set(headers
    manifold.h
    system.h
    system_extstep.h
    particle2d.h
    unicycle.h
    car.h
    rccar.h
    gunicycle.h
    body2d.h
    body3d.h
    uuv.h
    joint.h
    qrotor.h
    hrotor.h
    heli.h
    kinbody2d.h
    wrench.h
    mbs.h
    chain.h
    chain1.h
    airm.h
    airbot.h
    airbase.h
    snake.h
    autodiff.h
    dem.h
    force.h
    body2dforce.h
    arm.h
    urdf_parser.h
    ins.h
    point3d.h
    imu.h
)

add_subdirectory (manifolds)
add_subdirectory (costs)
add_subdirectory (parametrization)
add_subdirectory (sensors)

append_subdir_files(sources sources manifolds)
append_subdir_files(sources sources costs)
append_subdir_files(sources sources parametrization)
append_subdir_files(sources sources sensors)

append_subdir_files(headers headers costs)
append_subdir_files(headers headers manifolds)
append_subdir_files(headers headers parametrization)
append_subdir_files(headers headers sensors)
