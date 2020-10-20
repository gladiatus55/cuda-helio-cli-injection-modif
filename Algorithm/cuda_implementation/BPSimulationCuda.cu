#include <stdio.h>
#include <iostream>
#include <chrono>
#include <ctime>
#include <curand.h>
#include <curand_kernel.h>
#include <string>

#include "../../Constants/CosmicConstants.cuh"
#include "../../Constants/ParamsCarrier.cuh"
#include "BPDefines.cuh"
#include "cudaErrorCheck.cuh"
#include "../../libs/FMT/format.h"
#include "../../libs/FMT/printf.h"

extern "C" void runBPMethod(simulationInputBP *simulation);
extern "C" int mkdirAndchdir(std::string directoryName);

__device__ int countAtomicBP;

__device__ float getTkininjBP(unsigned long long state) {
	unsigned long long  ownState = state;
	int modulo;
	if (ownState > (150 * 10000)) {
		ownState -= (__double2ll_rd(ownState / (150 * 10000)) * (150 * 10000));
	}
	if (ownState >= 10000) {
		modulo = __double2int_rd(ownState / 10000);
	}
	else {
		modulo = 0;
	}
	return ((modulo)+((ownState - (modulo * 10000) + 1) / 10000.0f));
}

__global__ void nullCountBP() {
	countAtomicBP= 0;
}

__global__ void curandInitializationBP(curandState_t *state) {
	int execID = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(clock(), execID, 0, &state[execID]);
}

__global__ void wCalcBP(float *Tkininj, float *pinj, int padding) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	float Tkin = getTkininjBP(BLOCK_SIZE_BP * THREAD_SIZE_BP * padding + id);
	float Rig = sqrtf(Tkin*(Tkin + (2 * T0)));
	float p = Rig * 1e9 * q / c;
	pinj[id] = p;
	Tkininj[id] = Tkin; 
}

__global__ void trajectorySimulationBP(float *pinj, trajectoryHistoryBP* history, int padding, curandState* state) {
	extern __shared__ int sharedMemory[];
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int idx = threadIdx.x;
	float r = 1.0f;
	float p = pinj[id];
	float beta, Rig, dr, pp;
	float Tkin = getTkininjBP(BLOCK_SIZE_BP * THREAD_SIZE_BP * padding + id);
	float2* generated = (float2*) sharedMemory;
	curandState* cuState = (curandState *)(&generated[THREAD_SIZE_BP]);
	cuState[idx] = state[blockIdx.x * blockDim.x + threadIdx.x];
	int count;
	bool generate = true;
	for (; r < 100.0002f;) {
		beta = sqrtf(Tkin*(Tkin + T0 + T0)) / (Tkin + T0);
		Rig = (p * c / q)/1e9;
		pp = p;
		p -= (2.0f*V*pp*dt / (3.0f*r));
		if (generate) {
			generated[idx] = curand_box_muller(&cuState[idx]);
			dr = (V + (2.0f*K0*beta*Rig / r))*dt + (generated[idx].x * sqrtf(2.0f*K0*beta*Rig *dt));
			r += dr;
			generate = false;
		}
		else {
			dr = (V + (2.0f*K0*beta*Rig / r))*dt + (generated[idx].y * sqrtf(2.0f*K0*beta*Rig *dt));
			r += dr;
			generate = true;
		}
		Rig = p * c / q;
		Tkin = (sqrtf((T0*T0*q*q*1e9f*1e9f) + (q*q*Rig *Rig)) - (T0*q*1e9f)) / (q*1e9f);
		Rig = Rig / 1e9;
		beta = sqrtf(Tkin*(Tkin + T0 + T0)) / (Tkin + T0);
		if (beta>0.01f&&Tkin<200.0f) {
			if ((r > 100.0f) && ((r - dr) < 100.0f)) {
				count = atomicAdd(&countAtomicBP, 1);
				double newW = (m0_double*m0_double*c_double*c_double*c_double*c_double) + (p*p*c_double*c_double);
				newW = (pow(newW, -1.85) / p) / 1e45;
				history[count].setValues(Tkin, r, newW, id);
				break; 
			}
		}
		else if (beta<0.01f) {
			break;
		}
		if (r<0.3f) {
			r -= dr;
			p = pp;
		}
	}
	state[id] = cuState[idx];
}

void setConstantsBP(ParamsCarrier *singleTone) {
	float newDT = singleTone->getFloat("dt", -1.0f);
	if (newDT != -1.0f) {
		cudaMemcpyToSymbol(dt, &newDT, sizeof(newDT));
	}
	float newK = singleTone->getFloat("K0", -1.0f);
	if (newK != -1.0f) {
		cudaMemcpyToSymbol(K0, &newK, sizeof(newK));
	}
	float newV = singleTone->getFloat("V", 1.0f) * (-1.0f);
	if (newK != -1.0f) {
		cudaMemcpyToSymbol(V, &newV, sizeof(newV));
	}else{
		newV = 2.66667e-6 * (-1.0f); 		
		cudaMemcpyToSymbol(V, &newV, sizeof(newV));
	}
}

void runBPMethod(simulationInputBP *simulation) {
	fmt::MemoryWriter writer;
	int counter;
	ParamsCarrier *singleTone;
	singleTone = simulation->singleTone;
 	#if defined(_WIN32)
    		printf("Windows is not supported!"); 
		return; 
     	#else 
		mode_t process_mask = umask(0);		
		umask(process_mask);
		if(mkdirAndchdir("results") != 1){
			return; 
		} 
		if(mkdirAndchdir("BP") != 1){
			return; 
		}
		std::string destination = singleTone->getString("destination","");
		if(destination.length() == 0) {
			time_t clk = time(NULL);
			char *dirName = ctime(&clk); 
			for(int i = 0 ; i < strlen(dirName); i++){
				if(dirName[i] == ' '){
					dirName[i] = '_'; 
				}else if(i == strlen(dirName) - 1){
					dirName[i] = 0; 
				}
			}
			std::string parameters; 
			parameters += dirName; 
			parameters += "_dt=" + std::to_string(singleTone->getFloat("dt", 5.0f)) + "K0=" + 
				singleTone->getString("K0input", "0.000222") + "V=" + singleTone->getString("Vinput", "400"); 
			if(mkdirAndchdir(parameters.c_str()) != 1){
				return; 
			}
		} else {
			if(mkdirAndchdir(destination.c_str()) != 1){
				return; 
			}
		}
	#endif
	FILE  *bb = fopen("log.dat", "w");
	curandInitializationBP << <simulation->blockSize, simulation->threadSize >> >(simulation->state);
	gpuErrchk(cudaDeviceSynchronize());
	int iterations = singleTone->getInt("bilions", 1);
	setConstantsBP(singleTone);
	if (simulation->threadSize == 1024) {
		cudaFuncSetAttribute(trajectorySimulationBP, cudaFuncAttributeMaxDynamicSharedMemorySize, 65536);
	}
	for (int k = 0; k < iterations; ++k) {
		for (int i = 0; i < 1; ++i) {
			nullCountBP << <1, 1 >> >();
			gpuErrchk(cudaDeviceSynchronize());
			wCalcBP << <simulation->blockSize, simulation->threadSize >> >(simulation->Tkininj, simulation->pinj, i);
			gpuErrchk(cudaDeviceSynchronize());
			trajectorySimulationBP << <simulation->blockSize, simulation->threadSize, simulation->threadSize * sizeof(curandState_t) + 
				simulation->threadSize * sizeof(float2) >> >(simulation->pinj, simulation->history, i, simulation->state);
			gpuErrchk(cudaDeviceSynchronize());
			cudaMemcpyFromSymbol(&counter, countAtomicBP, sizeof(int), 0, cudaMemcpyDeviceToHost);
			if (counter != 0) {
				gpuErrchk(cudaMemcpy(simulation->local_history, simulation->history, counter*sizeof(trajectoryHistoryBP), cudaMemcpyDeviceToHost));
				for (int j = 0; j < counter; ++j) {
					fmt::fprintf(bb, "%g %g %g %g\n", simulation->local_history[j].Tkin, simulation->Tkininj[simulation->local_history[j].id] , simulation->local_history[j].r, simulation->local_history[j].w);
				}
			}
		}
	}
	fclose(bb);
}