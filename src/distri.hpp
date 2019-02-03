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

#include "typedef.hpp"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <memory>

using namespace std;

//------------------------------------------------------------------------------
#ifndef __DISTRI_HPP
#define __DISTRI_HPP
//------------------------------------------------------------------------------



struct border{

};

class Distri {
public:
  /// Constructs a grid based on a geometry
  Distri(const Geometry *geom);

  /// Constructs a grid based on a geometry with an offset
  // @param geom   Geometry information
  // @param offset distance of staggered grid point to cell's anchor point;
  //               (anchor point = lower left corner)
  Distri(const Geometry *geom, const multi_real_t &offset);

  /// Deletes the grid
  ~Distri();

  ///	Initializes the grid with a value
  void Initialize(const real_t &value);

  /// Write access to the grid cell at position [it]
  real_t &Cell(const Iterator &it, index_t vel);
  /// Read access to the grid cell at position [it]
  const real_t &Cell(const Iterator &it, index_t vel) const;
  const real_t &t(index_t vel) const;
  const real_t &ex(index_t vel) const;
  const real_t &ey(index_t vel) const;
  const real_t &rho(const Iterator &it) const;
  real_t &rho(const Iterator &it);
  real_t sum_vel(const Iterator &it);
  real_t sum_c_vel_x(const Iterator &it);
  real_t sum_c_vel_y(const Iterator &it);


private:
  real_t *_f1;
  real_t *_f2;
  real_t *_f3;
  real_t *_f4;
  real_t *_f5;
  real_t *_f6;
  real_t *_f7;
  real_t *_f8;
  real_t *_f0;
  real_t *_ex;
  real_t *_ey;
  real_t *_t;
  real_t *_rho;
  const Geometry *_geom;
};
//------------------------------------------------------------------------------
#endif // __GRID_HPP
