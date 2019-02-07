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

#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "parameter.hpp"
#include <stdio.h>
#include <string.h>



Geometry::~Geometry() {
    if (_cell)
        delete[] _cell;
}
//------------------------------------------------------------------------------
void Geometry::Load(const char *file) {
   FILE *handle = fopen(file, "r");
   real_t inval[2];
   real_t inval_vel[8];
   index_t inval_ind[4];
   char name[20000000];
   while (!feof(handle)) {
     if (!fscanf(handle, "%s =", name))
       continue;
     if (strcmp(name, "size") == 0) {
       if (fscanf(handle, " %f %f\n", &inval[0], &inval[1])) {
         _size[0] = inval[0];
         _size[1] = inval[1];
         printf("size: %d, %d \n", _size[0], _size[1]);
       }
       continue;
     }
     if (strcmp(name, "length") == 0) {
       if (fscanf(handle, " %f %f\n", &inval[0], &inval[1])) {
         _length[0] = inval[0];
         _length[1] = inval[1];
       }
       continue;
     }
     if (strcmp(name, "velocity") == 0) {
       if (fscanf(handle, " %f %f %f %f %f %f %f %f\n", &inval_vel[0], &inval_vel[1], &inval_vel[2], &inval_vel[3], &inval_vel[4], &inval_vel[5], &inval_vel[6], &inval_vel[7])) {
         _velocity[0] = inval_vel[0];
         _velocity[1] = inval_vel[1];
         for (int i = 0; i < 8; ++i){
        	 _bound_vel[i] = inval_vel[i];
         }
       }
       continue;
     }
     if (strcmp(name, "boundaries") == 0) {
       if (fscanf(handle, " %u %u %u %u\n", &inval_ind[0], &inval_ind[1], &inval_ind[2], &inval_ind[3])) {
         for (int i = 0; i < 4; ++i){
        	 _bound_true[i] = inval_ind[i];
         }
       }
       continue;
     }
     if (strcmp(name, "pressure") == 0) {
       if (fscanf(handle, " %f\n", &inval[0]))
         _pressure = inval[0];
       continue;
     }
     if (strcmp(name, "geometry") == 0) {
       if (!fscanf(handle, "%s\n", name))
         continue;
       if (strcmp(name, "free") == 0) {
         if (_cell)
           delete[] _cell;
         _cell = new Cell_t[_size[0] * _size[1]];
         bool parabolic = false;
         // Read stuff from file
         for (int y = _size[1]; y-- > 0;) {
           if (feof(handle)) {
             delete[] _cell;
             _cell = NULL;
             break;
           }
           for (int x = 0; x < _size[0]; ++x) {
             _cell[x + y * _size[0]].fluid = cellNone;
             _cell[x + y * _size[0]].factor = 1.0;
             switch (getc(handle)) {
             case '#':
               _cell[x + y * _size[0]].type = typeSolid;
               break;
             case 'I':
               _cell[x + y * _size[0]].type = typeIn;
               break;
             case 'O':
               _cell[x + y * _size[0]].type = typeOut;
               break;
             case '|':
               _cell[x + y * _size[0]].type = typeSlipV;
               break;
             case 'T':
               _cell[x + y * _size[0]].type = typeHot;
               break;
             case 'C':
               _cell[x + y * _size[0]].type = typeCold;
               break;
             case 'W':
               _cell[x + y * _size[0]].type = typeInsul;
               break;
             case '-':
               _cell[x + y * _size[0]].type = typeSlipH;
               break;
             case 'H':
               _cell[x + y * _size[0]].type = typeInH;
               break;
             case 'V':
               _cell[x + y * _size[0]].type = typeInV;
               break;
             default:
               if (x == 0 || x == _size[0] - 1 || y == 0 || y == _size[1] - 1)
                 _cell[x + y * _size[0]].type = typeSolid;
               else
                 _cell[x + y * _size[0]].type = typeFluid;
               break;
             };
           }
           if (!fscanf(handle, "\n"))
             continue;
         }
         if (!_cell)
           break;

         _size[0] -= 2;
         _size[1] -= 2;
       }
     }
   }
   fclose(handle);
   _h[0] = _length[0] / _size[0];
   _h[1] = _length[1] / _size[1];

   _blength = _length;

   _boffset = 0;



   _size[0] += 2;
   _size[1] += 2;

    _bsize = _size;
 }
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/*Geometry::Geometry(const Communicator *comm) : _comm(comm) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;

  _bsize[0] = _size[0] / _comm->ThreadDim()[0] + 2;
  _bsize[1] = _size[1] / _comm->ThreadDim()[1] + 2;
  if (_comm->ThreadIdx()[0] == _comm->ThreadDim()[0] - 1)
    _bsize[0] += _size[0] % _comm->ThreadDim()[0];
  if (_comm->ThreadIdx()[1] == _comm->ThreadDim()[1] - 1)
    _bsize[1] += _size[1] % _comm->ThreadDim()[0];

  _blength[0] = _h[0] * (_bsize[0] - 2);
  _blength[1] = _h[1] * (_bsize[1] - 2);

  // create boundary halo
  _size[0] += 2;
  _size[1] += 2;

  _cell = NULL;
  _boffset = 0;
}*/
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
Geometry::Geometry() {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;

  // create boundary halo
  _size[0] += 2;
  _size[1] += 2;

  _bsize = _size;
  _blength = _length;

  _cell = NULL;
  _boffset = 0;
}


//------------------------------------------------------------------------------
void Geometry::Update_Distr(Distri *f, const real_t &value,
                                     const Iterator &it) const {
  switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].fluid) {
  // cellW means: left neighbor is a fluid element
  case !cellW:
    break;
  case cellNW:
    break;
  case cellN:
    break;
  case cellNE:
    break;
  case cellE:
    break;
  case cellSE:
    break;
  case cellS:
    break;
  case cellSW:
    break;
  };
}
//------------------------------------------------------------------------------

const multi_index_t &Geometry::Size() const { return _bsize; }
//------------------------------------------------------------------------------
const multi_index_t &Geometry::TotalSize() const { return _size; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::Length() const { return _blength; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::TotalLength() const { return _length; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::Mesh() const { return _h; }
//------------------------------------------------------------------------------
void Geometry::Update_Fluid(Grid *u) const {
  if (_cell) {
    Iterator it(this);
    for (it.First(); it.Valid(); it.Next()) {
      // Update_Distr(u,0.0,it);
    }
  }
}


//------------------------------------------------------------------------------
