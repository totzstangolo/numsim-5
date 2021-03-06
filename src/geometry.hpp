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
//#include "communicator.hpp"
#include "iterator.hpp"
#include "parameter.hpp"
//------------------------------------------------------------------------------
#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP
//------------------------------------------------------------------------------
/// Typedef for cell types
typedef enum {
  typeFluid, // Standard fluid cell
  typeSolid, // Simple wall, no slip
  typeHot,
  typeCold,
  typeInsul,
  typeIn,    // Simple inflow (forced Velocity)
  typeOut,   // Outflow
  typeSlipH, // Horizontal slip boundary
  typeSlipV, // Vertical slip boundary
  typeInH,   // Horizontal inflow (parabolic)
  typeInV    // Vertical inflow (parabolic)
} CellType_t;

typedef enum {
  cellNone = 0, // Cell is not surrounded by fluid, or is fluid cell
  cellN = 1,    // Cell N is fluid
  cellW = 2,    // Cell W is fluid
  cellNW = 3,   // Cells N & W are fluid
  cellS = 4,    // Cell S is fluid
  cellSW = 6,   // Cells S & W are fluid
  cellE = 8,    // Cell E is fluid
  cellNE = 9,   // Cells N & E are fluid
  cellSE = 12   // Cells S & E are fluid
} CellBoundary_t;
/// Typedef for cells
typedef struct {
  CellType_t type;      // Cell type
  CellBoundary_t fluid; // Surrounding fluid cells
//  CellTemp_t temp; 		// Temperature Boundary Type
  real_t factor;        // Scale factor for inital values
} Cell_t;
//------------------------------------------------------------------------------
class Geometry {
public:
  /// Constructs a default geometry:
  // driven cavity with 128 x 128 grid, no-slip boundary conditions
  // as shown below
  //
  //       u=1, v=0
  //    -------------
  //    |           |
  // u=0|           |u=0
  // v=0|           |v=0
  //    |           |
  //    |           |
  //    -------------
  //       u=0, v=0
  Geometry();
  ~Geometry();

  /// Loads a geometry from a file
  void Load(const char *file);

  /// Returns the number of cells in each dimension
  const multi_index_t &Size() const;
  /// Returns the total number of cells in each dimension
  const multi_index_t &TotalSize() const;
  /// Returns the length of the domain
  const multi_real_t &Length() const;
  /// Returns the total length of the domain
  const multi_real_t &TotalLength() const;
  /// Returns the meshwidth
  const multi_real_t &Mesh() const;


  /// Returns the celltype of cell
  const CellType_t &get_cellType(index_t index) const {return _cell[index].type;};


  /// Updates the velocity field u
  void Update_Fluid(Grid *u) const;
  /// Updates the velocity field v
  void Update_V(Grid *v) const;
  /// Updates the pressure field p
  void Update_P(Grid *p) const;
  /// Updates the pressure field t
  void Update_T(Grid *t, real_t hot, real_t cold) const;

  void Update_Distr(Distri *f, const real_t &value,
                                       const Iterator &it) const;

  multi_real_t GetInflowVelo() const{
      return _velocity;
  }

  index_t* GetBoundStat() const{
        return _bound_true;
    }

  real_t* GetBoundVel() const{
      return _bound_vel;
  }

  real_t *_bound_vel = new real_t[8];
  index_t *_bound_true = new index_t[4];

private:
  const index_t _comm = 0;

  Cell_t *_cell;

  multi_index_t _size;
  multi_index_t _bsize;
  multi_real_t _length;
  multi_real_t _blength;
  index_t _boffset;
  multi_real_t _h;

  multi_real_t _velocity;
  real_t _pressure;
  real_t _T;
  multi_real_t _T_wall;

  void UpdateCellDirichlet_U(Grid *u, const real_t &value,
                             const Iterator &it) const;
  void UpdateCellDirichlet_V(Grid *v, const real_t &value,
                             const Iterator &it) const;
  void UpdateCellDirichlet_T(Grid *t, const real_t &value,
                             const Iterator &it) const;
  void UpdateCellNeumann(Grid *grid, const Iterator &it) const;
  void UpdateCellNeumann_P(Grid *grid, const Iterator &it) const;
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP
