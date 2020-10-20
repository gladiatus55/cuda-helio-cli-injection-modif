#ifndef COSMIC_CONSTANTS_H
#define COSMIC_CONSTANTS_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

static __device__ __constant__ float V = 2.66667e-6;    
__device__ __constant__ float dt = 5.0;
__device__ __constant__ float K0 = 0.000222;   
static __device__ __constant__ float m0 = 1.67261e-27;
static __device__ __constant__ float q = 1.60219e-19;
static __device__ __constant__ double m0_double = 1.67261e-27;
static __device__ __constant__ double q_double = 1.60219e-19;
static __device__ __constant__ double q_double_pow = 1.60219e-19 * 1.60219e-19;
static __device__ __constant__ double c_double = 2.99793e8;
static __device__ __constant__ float c = 2.99793e8;
static __device__ __constant__ float T0 = 1.67261e-27 * 2.99793e8 * 2.99793e8 / (1.60219e-19*1e9);
static __device__ __constant__ double T0_double = 1.67261e-27 * 2.99793e8 * 2.99793e8 / (1.60219e-19*1e9);
static __device__ __constant__ float T0w = 1.67261e-27 * 2.99793e8 * 2.99793e8;
static __device__ __constant__ float kparKper = 0.01f;
static __device__ __constant__ float omega2 = 7.25445387e-12; 
static __device__ __constant__ float Pi = 3.141592654;

void setDT(float newDT); 
void setK0(float newK0);

#endif 
