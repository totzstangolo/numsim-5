#include <assert.h>
#include <cuda.h>
#include <stdio.h>
#include "typedef.hpp"
#include "compute.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <unistd.h>
#include <vector>
#include <stdlib.h>
#include <memory>

real_t *d_w = nullptr;
real_t *d_index_x = nullptr;
real_t *d_index_y = nullptr;
real_t *d_bound_vel = nullptr;
index_t *d_bound_stat = nullptr;
real_t *d_boundary = nullptr;
real_t *d_rho = nullptr;
real_t *d_u = nullptr;
real_t *d_v = nullptr;
real_t *d_f = nullptr;
int *d_inv = nullptr;
int *d_x_delta = nullptr;
int *d_y_delta = nullptr;
//Data *d_grid = nullptr;

void InitGpu(Data *grid, index_t n){
	// Offload Data
	CUDA_CALL(cudaMalloc(&d_w, sizeof(real_t) * 9));
	CUDA_CALL(cudaMemcpy(d_w, &(grid->w[0]),
			sizeof(real_t) * 9, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_index_x, sizeof(real_t) * 9));
	CUDA_CALL(cudaMemcpy(d_index_x, &(grid->index_x[0]),
			sizeof(real_t) * 9, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_index_y, sizeof(real_t) * 9));
	CUDA_CALL(cudaMemcpy(d_index_y, &(grid->index_y[0]),
			sizeof(real_t) * 9, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_bound_vel, sizeof(real_t) * 8));
	CUDA_CALL(cudaMemcpy(d_bound_vel, &(grid->bound_vel[0]),
			sizeof(real_t) * 8, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_bound_stat, sizeof(index_t) * 4));
	CUDA_CALL(cudaMemcpy(d_bound_stat, &(grid->bound_stat[0]),
			sizeof(index_t) * 4, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_boundary, sizeof(real_t) * n));
	CUDA_CALL(cudaMemcpy(d_boundary, &(grid->boundary[0]),
			sizeof(real_t) * n, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_rho, sizeof(real_t) * n));
	CUDA_CALL(cudaMemcpy(d_rho, &(grid->rho[0]),
			sizeof(real_t) * n, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_u, sizeof(real_t) * n));
	CUDA_CALL(cudaMemcpy(d_u, &(grid->u[0]),
			sizeof(real_t) * n, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_v, sizeof(real_t) * n));
	CUDA_CALL(cudaMemcpy(d_v, &(grid->v[0]),
			sizeof(real_t) * n, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_f, sizeof(real_t) * n * 18));
	CUDA_CALL(cudaMemcpy(d_f, &(grid->f[0][0]),
			sizeof(real_t) * n * 18, cudaMemcpyHostToDevice));

	int inv[9] = {0, 2, 1, 4, 3, 6, 5, 8, 7};
	int x_delta[9] = {0, 1, -1, 0, 0, 1, -1, 1, -1};
	int y_delta[9] = {0, 0, 0, 1, -1, 1, -1, -1, 1};

	CUDA_CALL(cudaMalloc(&d_inv, sizeof(int) * 9));
	CUDA_CALL(cudaMemcpy(d_inv, &inv[0],
			sizeof(int) * 9, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_x_delta, sizeof(int) * 9));
	CUDA_CALL(cudaMemcpy(d_x_delta, &x_delta[0],
			sizeof(int) * 9, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_y_delta, sizeof(int) * 9));
	CUDA_CALL(cudaMemcpy(d_y_delta, &y_delta[0],
			sizeof(int) * 9, cudaMemcpyHostToDevice));


	/*d_grid = new Data(n);


	CUDA_CALL(cudaMalloc(&(d_grid->w), sizeof(real_t) * 9));
	CUDA_CALL(cudaMemcpy(d_grid->w, &(grid->w[0]),
			sizeof(real_t) * 9, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_grid->index_x, sizeof(real_t) * 9));
	CUDA_CALL(cudaMemcpy(d_grid->index_x, &(grid->index_x[0]),
			sizeof(real_t) * 9, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_grid->index_y, sizeof(real_t) * 9));
	CUDA_CALL(cudaMemcpy(d_grid->index_y, &(grid->index_y[0]),
			sizeof(real_t) * 9, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_grid->bound_vel, sizeof(real_t) * 8));
	CUDA_CALL(cudaMemcpy(d_grid->bound_vel, &(grid->bound_vel[0]),
			sizeof(real_t) * 8, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_grid->bound_stat, sizeof(index_t) * 4));
	CUDA_CALL(cudaMemcpy(d_grid->bound_stat, &(grid->bound_stat[0]),
			sizeof(index_t) * 4, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_grid->boundary, sizeof(real_t) * n));
	CUDA_CALL(cudaMemcpy(d_grid->boundary, &(grid->boundary[0]),
			sizeof(real_t) * n, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_grid->rho, sizeof(real_t) * n));
	CUDA_CALL(cudaMemcpy(d_grid->rho, &(grid->rho[0]),
			sizeof(real_t) * n, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_grid->u, sizeof(real_t) * n));
	CUDA_CALL(cudaMemcpy(d_grid->u, &(grid->u[0]),
			sizeof(real_t) * n, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_grid->v, sizeof(real_t) * n));
	CUDA_CALL(cudaMemcpy(d_grid->v, &(grid->v[0]),
			sizeof(real_t) * n, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_grid->f, sizeof(real_t) * n * 18));
	CUDA_CALL(cudaMemcpy(d_grid->f, &(grid->f[0]),
			sizeof(real_t) * n * 18, cudaMemcpyHostToDevice));*/


	printf("hi \n");
}

void KernelLaunch(index_t n, multi_index_t m, Data *grid, real_t omega){


	size_t threads_per_block = 256;
	size_t number_of_blocks = (n + threads_per_block - 1) / threads_per_block;

	size_t threads_small = 32;
	size_t blocks_small_x = (m[0] + threads_small - 1) / threads_small;
	size_t blocks_small_y = (m[1] + threads_small - 1) / threads_small;


	ForcesKernel<<<number_of_blocks, threads_per_block>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f, d_bound_vel);
	cudaDeviceSynchronize();

	if(grid->bound_stat[0]){
		BoundKernelEast<<<blocks_small_y, threads_small>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f, d_bound_vel);
	}
		cudaDeviceSynchronize();

	if(grid->bound_stat[1]){
		BoundKernelWest<<<blocks_small_y, threads_small>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f, d_bound_vel);
	}
	cudaDeviceSynchronize();

	if(grid->bound_stat[2]){
		BoundKernelNorth<<<blocks_small_x, threads_small>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f, d_bound_vel);
	}
	cudaDeviceSynchronize();

	if(grid->bound_stat[3]){
		BoundKernelSouth<<<blocks_small_x, threads_small>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f, d_bound_vel);
	}
	cudaDeviceSynchronize();

	CollisionKernel<<<number_of_blocks, threads_per_block>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f, d_bound_vel,
			omega, d_w, d_index_x, d_index_y);
	cudaDeviceSynchronize();

	StreamingKernel<<<(number_of_blocks), threads_per_block>>>(n, m[0], m[1], d_boundary, d_f, d_inv, d_x_delta, d_y_delta);
}

void CopyToCpu(index_t n, real_t * u_tmp, real_t * v_tmp, real_t * p_tmp){
	CUDA_CALL(cudaMemcpy(u_tmp, d_u, sizeof(real_t) * n,cudaMemcpyDeviceToHost));
	CUDA_CALL(cudaMemcpy(v_tmp, d_v, sizeof(real_t) * n,cudaMemcpyDeviceToHost));
	CUDA_CALL(cudaMemcpy(p_tmp, d_rho, sizeof(real_t) * n,cudaMemcpyDeviceToHost));
}

void FreeCuda(){
	CUDA_CALL(cudaFree(d_w));
	CUDA_CALL(cudaFree(d_index_x));
	CUDA_CALL(cudaFree(d_index_y));
	CUDA_CALL(cudaFree(d_bound_vel));
	CUDA_CALL(cudaFree(d_bound_stat));
	CUDA_CALL(cudaFree(d_boundary));
	CUDA_CALL(cudaFree(d_rho));
	CUDA_CALL(cudaFree(d_u));
	CUDA_CALL(cudaFree(d_v));
	CUDA_CALL(cudaFree(d_f));
	CUDA_CALL(cudaFree(d_inv));
	CUDA_CALL(cudaFree(d_x_delta));
	CUDA_CALL(cudaFree(d_y_delta));
}



__global__ void ForcesKernel(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel){
	index_t idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < n){

			if(!d_boundary[idx]){
				d_rho[idx] = 0;

				// computing rho, u and v
				for(index_t j = 0; j < 9; ++j){
					d_rho[idx] += d_f[idx + j * n];
				}
				d_u[idx] = (d_f[idx + n] - d_f[idx + 2*n] + d_f[idx + 5*n] - d_f[idx + 6*n] + d_f[idx + 7*n] - d_f[idx + 8*n]) / d_rho[idx];
				d_v[idx] = (d_f[idx + 3*n] - d_f[idx + 4*n] + d_f[idx +5*n] - d_f[idx + 6*n] - d_f[idx + 7*n] + d_f[idx + 8 *n]) / d_rho[idx];
			}
		}
}

__global__ void BoundKernelEast(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel){
	// inflow from east
	index_t idx = (1 + threadIdx.x + blockIdx.x * blockDim.x) * m0 - 1;
	if(idx < n){
		d_u[idx] = d_bound_vel[0];
		d_v[idx] = d_bound_vel[1];
		d_rho[idx] = 1.0/(1.0 +d_bound_vel[0]) *
			(  d_f[idx] +  d_f[idx + 3*n] + d_f[idx + 4*n]
			 + 2.0 * (d_f[idx +n] + d_f[idx +5*n] + d_f[idx+7*n]));
		//
		d_f[idx + 2*n] = d_f[idx + n] - 2/3 *  d_rho[idx] * d_bound_vel[0];
		d_f[idx + 6*n] = d_f[idx + 5*n] + 1/2*(d_f[idx+ 3*n] - d_f[idx +4*n])
								- 1/2*d_rho[idx]*d_bound_vel[1] - 1/6*d_rho[idx]*d_bound_vel[0];
		d_f[idx + 8*n] = d_f[idx + 7*n] + 1/2*(d_f[idx + 4*n] - d_f[idx + 3*n])
								+ 1/2*d_rho[idx]*d_bound_vel[1] - 1/6*d_rho[idx]*d_bound_vel[0];
	}
}

__global__ void BoundKernelWest(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel){
	// inflow from west
	index_t idx = (threadIdx.x + blockIdx.x * blockDim.x) * m0;
	if(idx < n){
		d_u[idx] = d_bound_vel[2];
		d_v[idx] = d_bound_vel[3];
		d_rho[idx] = 1.0/(1.0 +d_bound_vel[2]) *
			(  d_f[idx] +  d_f[idx + 3*n] + d_f[idx + 4*n]
			 + 2.0 * (d_f[idx +2*n] + d_f[idx +6*n] + d_f[idx+8*n]));
		//
		d_f[idx + 1*n] = d_f[idx + 2*n] + 2/3 *  d_rho[idx] * d_bound_vel[2];
		d_f[idx + 5*n] = d_f[idx + 6*n] + 1/2*(d_f[idx+ 4*n] - d_f[idx +3*n])
								- 1/2*d_rho[idx]*d_bound_vel[3] + 1/6*d_rho[idx]*d_bound_vel[2];
		d_f[idx + 7*n] = d_f[idx + 8*n] + 1/2*(d_f[idx + 3*n] - d_f[idx + 4*n])
								- 1/2*d_rho[idx]*d_bound_vel[3] - 1/6*d_rho[idx]*d_bound_vel[2];
	}
}

__global__ void BoundKernelNorth(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel){
	// inflow from north
	index_t idx = n - m0 + (threadIdx.x + blockIdx.x * blockDim.x);
	if(idx < n){
		d_u[idx] = d_bound_vel[4];
		d_v[idx] = d_bound_vel[5];
		d_rho[idx] = 1.0/(1.0 +d_bound_vel[5]) *
			(  d_f[idx] +  d_f[idx + 1*n] + d_f[idx + 2*n]
			 + 2.0 * (d_f[idx +3*n] + d_f[idx +5*n] + d_f[idx+8*n]));
		//
		d_f[idx + 4*n] = d_f[idx + 3*n] - 2/3 *  d_rho[idx] * d_bound_vel[5];
		d_f[idx + 6*n] = d_f[idx + 5*n] + 1/2*(d_f[idx+ 1*n] - d_f[idx +2*n])
								- 1/2*d_rho[idx]*d_bound_vel[4] - 1/6*d_rho[idx]*d_bound_vel[5];
		d_f[idx + 7*n] = d_f[idx + 8*n] + 1/2*(d_f[idx + 2*n] - d_f[idx + 1*n])
								+ 1/2*d_rho[idx]*d_bound_vel[4] - 1/6*d_rho[idx]*d_bound_vel[5];
	}
}

__global__ void BoundKernelSouth(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel){
	// inflow from south
	index_t idx = (threadIdx.x + blockIdx.x * blockDim.x);
	if(idx < m0){
		d_u[idx] = d_bound_vel[6];
		d_v[idx] = d_bound_vel[7];
		d_rho[idx] = 1.0/(1.0 +d_bound_vel[7]) *
			(  d_f[idx] +  d_f[idx + 1*n] + d_f[idx + 2*n]
			 + 2.0 * (d_f[idx +4*n] + d_f[idx +6*n] + d_f[idx+7*n]));
		//
		d_f[idx + 3*n] = d_f[idx + 4*n] + 2/3 *  d_rho[idx] * d_bound_vel[7];
		d_f[idx + 5*n] = d_f[idx + 6*n] + 1/2*(d_f[idx+ 2*n] - d_f[idx +1*n])
								+ 1/2*d_rho[idx]*d_bound_vel[6] + 1/6*d_rho[idx]*d_bound_vel[7];
		d_f[idx + 8*n] = d_f[idx + 7*n] + 1/2*(d_f[idx + 1*n] - d_f[idx + 2*n])
								- 1/2*d_rho[idx]*d_bound_vel[6] + 1/6*d_rho[idx]*d_bound_vel[7];
	}
}

__global__ void CollisionKernel(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t *d_bound_vel,
		real_t omega, real_t* d_w, real_t* d_index_x, real_t* d_index_y){
	index_t idx = (threadIdx.x + blockIdx.x * blockDim.x);
	real_t cu,u2,v2, eq = 0.0;
	if(idx <n){
		if(!d_boundary[idx]){
			//collision
			u2 = d_u[idx] * d_u[idx];
			v2 = d_v[idx] * d_v[idx];
			for(index_t j = 0; j < 9; ++j){
				cu = 3.0 * (d_index_x[j] * d_u[idx] + d_index_y[j] * d_v[idx]);
				eq = d_rho[idx] * d_w[j] * (1 + cu + 0.5 * (cu * cu) - 1.5 * (u2 + v2));
				d_f[idx + (9 + j)*n] = d_f[idx + j*n] + omega * (eq - d_f[idx + j*n]);
			}
		}
	}
}


__global__ void StreamingKernel(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t *d_f, int * d_inv, int *d_x_delta, int *d_y_delta){
	index_t idx = (threadIdx.x + blockIdx.x * blockDim.x) + m0;
	if(idx < n - m0 && threadIdx.x != 0 && threadIdx.x != m0 - 1){
		// Boundarys don't stream
		if(!d_boundary[idx]){
			d_f[idx] = d_f[idx + 9 * n];
			for(int j = 1; j < 9; ++j){
				//neighbor is boundary
				//debug = k - x_delta[j] - y_delta[j] * m;
				if(d_boundary[idx - d_x_delta[j] - d_y_delta[j] * m0]){
					// bounceback
					d_f[idx + j*n] = d_f[idx + (d_inv[j] + 9) * n];
				}else{
					//neighbor is no boundary (standard)
					d_f[idx + j * n] = d_f[idx - d_x_delta[j] - d_y_delta[j] * m0 + (9 + j) * n];
				}
			}
		}
	}
}


__global__ void DebugKernel(Data * d_grid) {
	if(threadIdx.x == 0){
		printf("Hi \n");
		printf("w: %f \n", d_grid->w[0]);
	}
}
