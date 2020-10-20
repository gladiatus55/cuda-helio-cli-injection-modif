#ifndef BP_DEFINES_H
#define BP_DEFINES_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <curand.h>
#include <curand_kernel.h>

#if (__CUDA_ARCH__ == 610)
#define BLOCK_SIZE_BP 32768
#define THREAD_SIZE_BP 512
#elif (__CUDA_ARCH__ == 750)
#define  BLOCK_SIZE_BP 16384
#define THREAD_SIZE_BP 1024
#else 
#define BLOCK_SIZE_BP 64
#define THREAD_SIZE_BP 64
#endif

struct trajectoryHistoryBP {
	float Tkin = -1.0f;
	float r = -1.0f;
	double w = -1.0f;
	int id = -1;

	__device__ void setValues(float newTkin, float newR, float newW, int newId) {
		Tkin = newTkin;
		r = newR;
		w = newW;
		id = newId;
	}
};

struct simulationInputBP {
	ParamsCarrier *singleTone; 
	curandState_t *state;
	trajectoryHistoryBP *history;
	trajectoryHistoryBP *local_history; 
	float *Tkininj; 
	float *pinj; 
	double *w;
	int blockSize;
	int threadSize;
};


#endif // !BP_DEFINES_H
