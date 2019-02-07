# numsim-5

Go into build directory

    cmake  ../ && make

Alternatively try to compile with

    nvcc -std=c++11 -L/usr/local/cuda/lib64 -lcudart $(pkg-config --cflags --libs sdl2) -D USE_VTK ../src/*.cpp ../src/*.cu -o numsim

Stay in working directory. To run the cases, execute:

	time ./numsim -geom ../geom/cav128.geom -param ../param/def128.param

Remember to remove the VTK directory after each case run.
