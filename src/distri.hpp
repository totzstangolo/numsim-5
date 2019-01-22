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
//------------------------------------------------------------------------------
#ifndef __DISTRI_HPP
#define __DISTRI_HPP
//------------------------------------------------------------------------------
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
  // const real_t &Cell(const Iterator &it, index_t vel) const;


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
  const Geometry *_geom;
};
//------------------------------------------------------------------------------
#endif // __GRID_HPP
