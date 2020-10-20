#include "CosmicConstants.cuh"

#include <stdio.h>

__global__ void mykernel() {
	printf("value is %f\n", K0);
}

void setDT(float newDT) {
	cudaMemcpyToSymbol(dt, &newDT, sizeof(newDT));
}

void setK0(float newK0) {
	mykernel << <1, 1 >> >();
	cudaMemcpyToSymbol(K0, &newK0, sizeof(newK0), 0, cudaMemcpyHostToDevice);
	mykernel << <1, 1 >> >();
}

