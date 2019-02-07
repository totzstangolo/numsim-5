# Copyright (C) 2015   Michael Lahnert
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# add possibility to add a debug visu
vars = Variables('custom.py')
vars.Add (BoolVariable('visu', 'Set to 1 for enabling debug visu', 0))
vars.Add (BoolVariable('vtk', 'Set to 1 for enabling VTK', 0))
vars.Add (BoolVariable('opt_dt', 'Set to 1 for enabling optimized dt value instead of the provided by param file, but wihout secured stability', 0))
vars.Add (BoolVariable('opt_omega', 'Set to 1 for enabling optimized omega value instead of the provided by param file, but wihout secured stability', 0))

env = Environment(variables=vars)

# do debug build?
debug = ARGUMENTS.get('debug', 1)

# set the compiler.
# For using clang in parallel you have to set all flags by hand or define a
# macro similar to mpic++

# parallel
# env.Replace(CXX='mpic++')

env.Replace(CXX='nvcc')

# serial
# env.Replace(CXX='g++')

# define some general compiler flags
env.Append(
           CXXFLAGS=[
                     "-Wall",
                     "-Wextra",
                     "-pedantic",
                     "-std=c++14",
                     "-Wno-unused-parameter",
                     "-L/usr/local/cuda/lib64",
                     "-lcudart",
                     # "$(pkg-config --cflags --libs sdl2)"
           ]
          )

# add flags for debug and release build
if debug == 0:
    env['CXXFLAGS'] += ["-O3"]
else:
    env['CXXFLAGS'] += [
                        "-g3",
                        "-O0",
                       ]

# call SConscript to actually build the project after setting up the environment
env.SConscript("./SConscript", exports='env', variant_dir='./build', duplicate=0)
