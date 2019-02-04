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
#include <cuda.h>
#include <cuda_runtime.h>

//------------------------------------------------------------------------------
#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP
//------------------------------------------------------------------------------

extern cudaError_t CU_err;
extern cudaError_t _err;


// Simple error checking function for CUDA actions
#define CUDA_CALL( call )               \
{                                       \
cudaError_t result = call;              \
if ( cudaSuccess != result )            \
    std::cerr << "CUDA error " << result << " in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString( result ) << " (" << #call << ")" << std::endl;  \
}


//------------------------------------------------------------------------------


struct Data{
	index_t sizeTest;
	real_t* w;
	real_t* index_x;
	real_t* index_y;
	real_t** f;
	real_t* boundary;
	real_t* rho;
	real_t* u;
	real_t* v;
	index_t* bound_stat;
	real_t* bound_vel;

	//constructor
	Data(const index_t size);

	//destructor
	~Data();

};

void InitGpu(Data *grid, index_t n);

void KernelLaunch(index_t n, multi_index_t m, Data *grid, real_t omega);

void CopyToCpu(index_t n, real_t * u_tmp, real_t * v_tmp, real_t * p_tmp);

void FreeCuda();

__global__ void ForcesKernel(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel);

__global__ void BoundKernelEast(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel);

__global__ void BoundKernelWest(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel);

__global__ void BoundKernelNorth(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel);

__global__ void BoundKernelSouth(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel);

__global__ void CollisionKernel(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel,
		real_t omega, real_t* d_w, real_t* d_index_x, real_t* d_index_y);

__global__ void StreamingKernel(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t *d_f, int * d_inv, int *d_x_delta, int *d_y_delta);



__global__ void DebugKernel(Data* grid);



class Compute {
public:
  /// Creates a compute instance with given geometry and parameter
  Compute(const Geometry *geom, const Parameter *param);
  /// Deletes all grids
  ~Compute();

  /// Execute one time step of the fluid simulation (with or without debug info)
  // @ param printInfo print information about current solver state (residual
  // etc.)
  void TimeStep(bool printInfo);

  void Init();

  /// Returns the simulated time in total
  const real_t &GetTime() const;

  /// Returns the pointer to U
  const Grid *GetU() const;
  /// Returns the pointer to V
  const Grid *GetV() const;
  /// Returns the pointer to P
  const Grid *GetP() const;
  /// Returns the pointer to RHS
  const Grid *GetT() const;
  /// Returns the pointer to T
  const Grid *GetRHS() const;

  /// Computes and returns the absolute velocity
  const Grid *GetVelocity();
  /// Computes and returns the vorticity
  const Grid *GetVorticity();
  /// Computes and returns the stream line values
  const Grid *GetStream();

  static void CudaFree();

  /*real_t* GetU_tmp() const{ return u_tmp;}
  real_t* GetV_tmp() const{ return v_tmp;}
  real_t* Getp_tmp() const{ return p_tmp;}*/

  // copy dummies
  real_t * u_tmp;
  real_t * v_tmp;
  real_t * p_tmp;



private:
  // testing
  index_t _iter_count;
  real_t _max_dt;




  // current timestep
  real_t _t;

  // donor-cell diffusion condition (p. 27)
  real_t _dtlimit;

  // limit for residual
  real_t _epslimit;

  // velocities
  Grid *_u;
  Grid *_v;

  // pressure and temperature
  Grid *_p;
  Grid *_T;

  // prel. vel
  Grid *_F;
  Grid *_G;

  // right-hand side
  Grid *_rhs;

  // container for interpolating whichever values
  Grid *_tmp;

  //Distri *_f;
  //Distri *_f_new;
  //Distri *_f_eq;
  Data *grid;

  Solver *_solver;

  const Geometry *_geom;
  const Parameter *_param;

  /// Compute the new velocites u,v
  void NewVelocities(const real_t &dt);
  /// Compute the temporary velocites F,G
  void MomentumEqu(const real_t &dt);
  /// Compute the modified temporary velocites ~F,~G
  void ModMomentumEqu(const real_t &dt);
  /// Compute the temporary velocites F,G
  void TempEqu(const real_t &dt);
  /// Compute the RHS of the poisson equation
  void RHS(const real_t &dt);
};
//------------------------------------------------------------------------------
#endif // __COMPUTE_HPP
