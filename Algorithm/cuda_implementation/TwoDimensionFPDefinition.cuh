#ifndef FLOAT_TWO_FP_DEFINES_H
#define FLOAT_TWO_FP_DEFINES_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <curand.h>
#include <curand_kernel.h>

#if (__CUDA_ARCH__ == 610)
#define BLOCK_SIZE 65536
#define THREAD_SIZE 512
#elif (__CUDA_ARCH__ == 750)
#define BLOCK_SIZE 32768
#define THREAD_SIZE 1024
#else 
#define BLOCK_SIZE 64
#define THREAD_SIZE 64
#endif

struct trajectoryHistoryTwoFP {
	float p = -1.0f;
	float r = -1.0f;
    float mu = -1.0f;
    float LBT = -1.0f;
	int id = -1;

	__device__ void setValues(float newP, float newR, int newId, float newMu, float newLBT) {
        p = newP;
        r = newR;
		id = newId;
        mu = newMu;
        LBT = newLBT;
	}
};

struct simulationInputTwoFP {
	ParamsCarrier *singleTone;
	curandState_t *state;
	trajectoryHistoryTwoFP *history;
	trajectoryHistoryTwoFP *local_history; 
    float *pinj;
    float *mu; 
	double *w;
	int blockSize;
	int threadSize;
};


#endif // !FLOAT_FW_DEFINES_H
