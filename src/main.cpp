/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//------------------------------------------------------------------------------
#include "typedef.hpp"
//#include "communicator.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "grid.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <unistd.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <memory>
//#include "compute.cu"


#ifdef USE_DEBUG_VISU
#include "visu.hpp"
#endif // USE_DEBUG_VISU
#ifdef USE_VTK
#include "vtk.hpp"
#include <sys/stat.h>
#endif // USE_VTK

#include "argvparser.hpp"

using namespace std;


int main(int argc, char **argv) {

  // Create parameter and geometry instances with default values
  //Communicator comm(&argc, &argv);
  Parameter param;
  printf("Hi \n");

  Geometry geom;
  printf("Hi \n");

  ARGVParser parser;
  parser.bind("-geom", [&geom](int ac, char **av) -> int {
	geom.Load(av[0]);
    return 1;
  });
  printf("Hi \n");

  parser.bind("-param", [&param](int ac, char **av) -> int {
    param.Load(av[0]);
    return 1;
  });
  printf("Hi \n");

  parser.exec(argc, argv);
  printf("Hi \n");

  // Create the fluid solver
  Compute comp(&geom, &param);

  printf("Hi \n");






#ifdef USE_VTK
    // check if folder "VTK" exists
    struct stat info;

    if (stat("VTK", &info) != 0) {
      system("mkdir VTK");
    }
#endif

// Create and initialize the visualization
#ifdef USE_DEBUG_VISU
  Renderer visu(geom.Length(), geom.Mesh());
  double ratio = geom.Length()[1]/geom.Length()[0];
  visu.Init(800/ comm.ThreadDim()[0], 800*ratio/ comm.ThreadDim()[1], comm.getRank() + 1);
#endif // USE_DEBUG_VISU

// #ifdef USE_VTK
//   // Create a VTK generator;
//   // use offset as the domain shift
//   multi_real_t offset;
//   offset[0] = comm.ThreadIdx()[0] * (geom.Mesh()[0] * (double)(geom.Size()[0] - 2));
//   offset[1] = comm.ThreadIdx()[1] * (geom.Mesh()[1] * (double)(geom.Size()[1] - 2));
//   VTK vtk(geom.Mesh(), geom.Size(), geom.TotalSize(), offset, comm.getRank(),
//           comm.getSize(), comm.ThreadDim());
// #endif

multi_index_t dim;
dim[0] = 1;
dim[1] = 1;
#ifdef USE_VTK
VTK vtk(geom.Mesh(), geom.Length(), geom.TotalLength(), 0,
      1, dim);
#endif

#ifdef USE_DEBUG_VISU
  const Grid *visugrid;

  visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU

  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend()) {
#ifdef USE_DEBUG_VISU
    // Render and check if window is closed
    switch (visu.Render(visugrid)) {
    case -1:
      return -1;
    case 0:
      visugrid = comp.GetVelocity();
      break;
    case 1:
      visugrid = comp.GetU();
      break;
    case 2:
      visugrid = comp.GetV();
      break;
    case 3:
      visugrid = comp.GetP();
      break;
    case 4:
      visugrid = comp.GetT();
    default:
      break;
    };
#endif // USE_DEBUG_VISU

#ifdef USE_VTK
    // Create VTK Files in the folder VTK
    vtk.Init("VTK/field");
    vtk.AddRank();
    vtk.AddCellField("Cell Velocity", comp.GetU(), comp.GetV());
    vtk.SwitchToPointData();
    vtk.AddPointField("Velocity", comp.GetU(), comp.GetV());
    vtk.AddPointScalar("Pressure", comp.GetP());
    //vtk.AddPointScalar("Temperature", comp.GetT());
    vtk.Finish();
#endif

    // suppress output on other nodes than rank 0
    bool printOnlyOnMaster = 1;
    printf("Hi \n");
    comp.TimeStep(printOnlyOnMaster);
  }
  comp.CudaFree;
  return 0;
}
