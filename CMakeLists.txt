#
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
#

cmake_minimum_required (VERSION 3.9)

project (NUMSIM LANGUAGES CXX CUDA)
set(CXX nvcc)

# SET ( CMAKE_CXX_FLAGS "-D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES" CACHE STRING "compile flags" FORCE)

add_definitions(-DUSE_VTK)
# find mpi
#find_package (MPI REQUIRED)
#INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${MPI_INCLUDE_PATH})




add_executable(numsim
#  ${CMAKE_CURRENT_SOURCE_DIR}/src/communicator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grid.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameter.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/typedef.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/zeitgeist.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/src/communicator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.cu
)

#  target_link_libraries (numsim ${MPI_LIBRARIES})
#  set_target_properties(numsim PROPERTIES
#    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#  set_target_properties (numsim PROPERTIES
#    LINK_FLAGS "${MPI_LINK_FLAGS}")



target_compile_features(numsim PUBLIC cxx_std_14)

set_target_properties( numsim
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON
)
