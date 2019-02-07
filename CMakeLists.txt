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
#set(CXX nvcc)

option (VTK_VISU
  "Enable VTK visualization" ON)

if (VTK_VISU)
  add_definitions(-DUSE_VTK)
endif (VTK_VISU)

add_executable(numsim
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
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.cu
)




target_compile_features(numsim PUBLIC cxx_std_14)

set_target_properties( numsim
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON
)





