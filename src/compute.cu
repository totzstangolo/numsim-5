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
#include <stdio.h>
#include <stdlib.h>
#include <memory>

//Data *d_grid = nullptr;

/*void InitGpu(Data grid, index_t n){
	// init grid
	CUDA_CALL(cudaMalloc(&d_grid, sizeof(Data)));
	CUDA_CALL(cudaMemcpy(d_grid, &grid,
	                         sizeof(Data), cudaMemcpyHostToDevice));

	//copy data
	CUDA_CALL(cudaMalloc(&(d_grid->w), sizeof(real_t) * 9));
	CUDA_CALL(cudaMemcpy(d_grid->w, &grid.w,
		                         sizeof(real_t) * 9, cudaMemcpyHostToDevice));

	//debugkernel<<<1, 32>>>(d_grid);


	printf("hi");
}*/



/*__global__ void debugkernel(Data *d_grid) {
	if(threadIdx.x == 1)
	printf("w %f %f %f", d_grid->w[0], d_grid->w[1], d_grid->w[2]);
}*/
