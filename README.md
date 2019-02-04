# numsim-5
Compile with

    scons -c && scons

Stay in working directory. To run the cases, execute:

    ./build/NumSim -geom ./geom/cavity.geom -param ./param/default.param

    ./build/NumSim -geom ./geom/karman.geom -param ./param/karman.param

    ./build/NumSim -geom ./geom/step.geom -param ./param/default.param

    ./build/NumSim -geom ./geom/trap.geom -param ./param/default.param

Remember to remove the VTK directory after each case run.
