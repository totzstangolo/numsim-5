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
#include <cuda_runtime.h>


__constant__ index_t d_size[3];

__constant__ real_t d_w[9];
__constant__ real_t d_index_x[9];
__constant__ real_t d_index_y[9];
__constant__ real_t d_bound_vel[8];
//__constant__ index_t d_bound_stat[4]; not necessary
__constant__ int d_inv[9];
__constant__ int d_x_delta[9];
__constant__ int d_y_delta[9];

real_t *d_boundary = nullptr;
real_t *d_rho = nullptr;
real_t *d_u = nullptr;
real_t *d_v = nullptr;
real_t *d_f = nullptr;
cudaStream_t streams[4];





void InitGpu(Data *grid, index_t n, index_t m0, index_t m1){

	int inv[9] = {0, 2, 1, 4, 3, 6, 5, 8, 7};
	int x_delta[9] = {0, 1, -1, 0, 0, 1, -1, 1, -1};
	int y_delta[9] = {0, 0, 0, 1, -1, 1, -1, -1, 1};

	index_t sizes[3];
	sizes[0] = n;
	sizes[1] = m0;
	sizes[2] = m1;

	// Offload/Init constants
	// sizes
	CUDA_CALL(cudaMemcpyToSymbol(d_size,(void*)&sizes,sizeof(index_t) * 3,0, cudaMemcpyHostToDevice));

	// w
	CUDA_CALL(cudaMemcpyToSymbol(d_w,(void*)&(grid->w[0]),9*sizeof(real_t),0, cudaMemcpyHostToDevice));
	// index_x
	CUDA_CALL(cudaMemcpyToSymbol(d_index_x,(void*)&(grid->index_x[0]),9*sizeof(real_t),0, cudaMemcpyHostToDevice));
	// index_y
	CUDA_CALL(cudaMemcpyToSymbol(d_index_y,(void*)&(grid->index_y[0]),9*sizeof(real_t),0, cudaMemcpyHostToDevice));
	// bound_vel
	CUDA_CALL(cudaMemcpyToSymbol(d_bound_vel,(void*)&(grid->bound_vel[0]),8*sizeof(real_t),0, cudaMemcpyHostToDevice));
	// bound_stat not necessary
	// inv
	CUDA_CALL(cudaMemcpyToSymbol(d_inv,(void*)&inv[0],9*sizeof(int),0, cudaMemcpyHostToDevice));
	// x_delta
	CUDA_CALL(cudaMemcpyToSymbol(d_x_delta,(void*)&x_delta[0],9*sizeof(int),0, cudaMemcpyHostToDevice));
	// y_delta
	CUDA_CALL(cudaMemcpyToSymbol(d_y_delta,(void*)&y_delta[0],9*sizeof(int),0, cudaMemcpyHostToDevice));

	cudaStream_t streams[4];
	for(int i = 0; i< 4; ++i){
		cudaStreamCreate(&streams[i]);
	}

	printf("5\n");


	// Offload Grid
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

	printf("6\n");


	CUDA_CALL(cudaMalloc(&d_f, sizeof(real_t) * n * 18));
	CUDA_CALL(cudaMemcpy(d_f, &(grid->f[0]),
			sizeof(real_t) * n * 18, cudaMemcpyHostToDevice));

	printf("7\n");


}

void KernelLaunch(index_t n, multi_index_t m, Data *grid, real_t omega){


	//size_t THREADS_PER_BLOCK = 256;
	size_t number_of_blocks = (n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

	//size_t THREADS_SMALL = 32;
	size_t blocks_small_x = (m[0] + THREADS_SMALL - 1) / THREADS_SMALL;
	size_t blocks_small_y = (m[1] + THREADS_SMALL - 1) / THREADS_SMALL;


	ForcesKernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f);

	if(grid->bound_stat[0]){
		BoundKernelEast<<<blocks_small_y, THREADS_SMALL, 0, streams[0]>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f);
	}


	if(grid->bound_stat[1]){
		BoundKernelWest<<<blocks_small_y, THREADS_SMALL, 0, streams[1]>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f);
	}

	if(grid->bound_stat[2]){
		BoundKernelNorth<<<blocks_small_x, THREADS_SMALL, 0, streams[2]>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f);
	}

	if(grid->bound_stat[3]){
		BoundKernelSouth<<<blocks_small_x, THREADS_SMALL, 0, streams[3]>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f);
	}
	cudaDeviceSynchronize();


	CollisionKernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(n, m[0], m[1], d_boundary, d_rho, d_u, d_v, d_f, omega);

	size_t number_of_blocks_stream = (n + THREADS_PER_BLOCK - 3) / THREADS_PER_BLOCK;

	StreamingKernel<<< number_of_blocks_stream, THREADS_PER_BLOCK>>>(n, m[0], m[1], d_boundary, d_f);
}

void CopyToCpu(index_t n, real_t * u_tmp, real_t * v_tmp, real_t * p_tmp){
	CUDA_CALL(cudaMemcpy(u_tmp, d_u, sizeof(real_t) * n,cudaMemcpyDeviceToHost));
	CUDA_CALL(cudaMemcpy(v_tmp, d_v, sizeof(real_t) * n,cudaMemcpyDeviceToHost));
	CUDA_CALL(cudaMemcpy(p_tmp, d_rho, sizeof(real_t) * n,cudaMemcpyDeviceToHost));
}

void FreeCuda(){
	CUDA_CALL(cudaFree(d_boundary));
	CUDA_CALL(cudaFree(d_rho));
	CUDA_CALL(cudaFree(d_u));
	CUDA_CALL(cudaFree(d_v));
	CUDA_CALL(cudaFree(d_f));
}



__global__ void ForcesKernel(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f){
	index_t idx = threadIdx.x + blockIdx.x * blockDim.x;
	if(idx < n){

			if(!d_boundary[idx]){
				d_rho[idx] = 0;

				// computing rho, u and v
#pragma unroll(9)
				for(index_t j = 0; j < 9; ++j){
					d_rho[idx] += d_f[idx + j * n];
				}
				d_u[idx] = (d_f[idx + n] - d_f[idx + 2*n] + d_f[idx + 5*n] - d_f[idx + 6*n] + d_f[idx + 7*n] - d_f[idx + 8*n]) / d_rho[idx];
				d_v[idx] = (d_f[idx + 3*n] - d_f[idx + 4*n] + d_f[idx +5*n] - d_f[idx + 6*n] - d_f[idx + 7*n] + d_f[idx + 8 *n]) / d_rho[idx];
			}
		}
}

__global__ void BoundKernelEast(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f){
	// inflow from east
	index_t idx = (2 + threadIdx.x + blockIdx.x * blockDim.x) * m0 - 1;
	if(idx < n - m0){
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

__global__ void BoundKernelWest(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f){
	// inflow from west
	index_t idx = (1 + threadIdx.x + blockIdx.x * blockDim.x) * m0;
	if(idx < n - m0){
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

__global__ void BoundKernelNorth(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f){
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

__global__ void BoundKernelSouth(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f){
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

__global__ void CollisionKernel(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t * d_rho, real_t * d_u, real_t * d_v, real_t *d_f, real_t omega){
	index_t idx = (threadIdx.x + blockIdx.x * blockDim.x);
	real_t cu,u2,v2, eq = 0.0;
	if(idx <n){
		if(!d_boundary[idx]){
			//collision
			u2 = d_u[idx] * d_u[idx];
			v2 = d_v[idx] * d_v[idx];
#pragma unroll(9)
			for(index_t j = 0; j < 9; ++j){
				cu = 3.0 * (d_index_x[j] * d_u[idx] + d_index_y[j] * d_v[idx]);
				eq = d_rho[idx] * d_w[j] * (1 + cu + 0.5 * (cu * cu) - 1.5 * (u2 + v2));
				d_f[idx + (9 + j)*n] = d_f[idx + j*n] + omega * (eq - d_f[idx + j*n]);
			}
		}
	}
}


__global__ void StreamingKernel(index_t n, index_t m0, index_t m1, real_t * d_boundary, real_t *d_f){
	index_t k = (m0 + threadIdx.x + blockIdx.x * blockDim.x);
	index_t offset = k % m0;
	if(k < n - m0 && offset != 0 && offset !=m0-1){
		index_t idx = k;

		// Boundarys don't stream
		if(!d_boundary[idx]){
			d_f[idx] = d_f[idx + 9 * m0 * m1];
#pragma unroll(9)
			for(int j = 1; j < 9; ++j){
				//neighbor is boundary
				//debug = k - x_delta[j] - y_delta[j] * m;
				if(d_boundary[idx - d_x_delta[j] - d_y_delta[j] * m0]){
					// bounceback
					d_f[idx + j*m0 * m1] = d_f[idx + (d_inv[j] + 9) * m0 * m1];
				}else{
					//neighbor is no boundary (standard)
					d_f[idx + j * m0 * m1] = d_f[idx - d_x_delta[j] - d_y_delta[j] * m0 + (9 + j) * m0 * m1];
				}
			}
		}
	}
}


__global__ void DebugKernel() {
	if(threadIdx.x == 0){
		printf("inv2: %d w3: %f \n", d_inv[2], d_w[3]);
	}
}
