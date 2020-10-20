#include <stdio.h>
#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <curand.h>
#include <curand_kernel.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "../../Constants/CosmicConstants.cuh"
#include "../../Constants/ParamsCarrier.cuh"
#include "FloatFWDefines.cuh"
#include "cudaErrorCheck.cuh"
#include "../../libs/FMT/format.h"
#include "../../libs/FMT/printf.h"

extern "C" void runFWMethod(simulationInput *simulation);
extern "C" int mkdirAndchdir(std::string directoryName);

__device__ int countAtomic;

__device__ float getTkininj(unsigned long long state) {
	unsigned long long  ownState = state;
	int modulo;
	if (ownState > (101 * 10000)) {
		ownState -= (__double2ll_rd(ownState / (101 * 10000)) * (101 * 10000));
	}
	if (ownState >= 10000) {
		modulo = __double2int_rd(ownState / 10000);
	}
	else {
		modulo = 0;
	}
	return ((modulo)+((ownState - (modulo * 10000) + 1) / 10000.0f));
}

__global__ void nullCount() {
	countAtomic = 0;
}

__global__ void curandInitialization(curandState_t *state) {
	int execID = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(clock(), execID, 0, &state[execID]);
}

__global__ void wCalc(double *w, float *pinj, int padding, trajectoryHistory *history) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	float Tkinw = getTkininj(BLOCK_SIZE * THREAD_SIZE * padding + id)*1e9f*q;
	float Rig = sqrtf(Tkinw*(Tkinw + (2 * T0w))) / q;
	float p = Rig*q / c;
	double newW = (m0_double*m0_double*c_double*c_double*c_double*c_double) + (p*p*c_double*c_double);
	newW = (pow(newW, -1.85) / p) / 1e45;
	w[id] = newW;
	pinj[id] = p;
}

__global__ void trajectorySimulation(float *pinj, trajectoryHistory* history, int padding, curandState* state) {
	extern __shared__ int sharedMemory[];
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int idx = threadIdx.x;
	float r = 100.0001f;
	float p = pinj[id];
	float beta, sumac = 0.0f, Rig, dr, pp;
	float Tkin = getTkininj(BLOCK_SIZE * THREAD_SIZE * padding + id);
	float2* generated = (float2*) sharedMemory;
	curandState* cuState = (curandState *)(&generated[THREAD_SIZE]);
	cuState[idx] = state[id];
	int count;
	bool generate = true;
	for (; r < 100.0002f;) {
		beta = sqrtf(Tkin*(Tkin + T0 + T0)) / (Tkin + T0);
		Rig = sqrtf(Tkin*(Tkin + (2.0f * T0)));
		sumac += ((4.0f * V / (3.0f*r))*dt);
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
		Rig = p*c / q;
		Tkin = (sqrtf((T0*T0*q*q*1e9f*1e9f) + (q*q*Rig *Rig)) - (T0*q*1e9f)) / (q*1e9f);
		if (beta>0.01f&&Tkin<100.0f) {
			if ((r - 1.0f) / ((r - dr) - 1.0f) < 0.0f) {
				count = atomicAdd(&countAtomic, 1);
				history[count].setValues(sumac, r, p, id);
			}
		}
		else if (beta<0.01f) {
			break;
		}
		if (r<0.1f) {
			r -= dr;
			p = pp;
		}
	}
	state[id] = cuState[idx];
}

void setConstants(ParamsCarrier *singleTone) {
	float newDT = singleTone->getFloat("dt", -1.0f);
	if (newDT != -1.0f) {
		cudaMemcpyToSymbol(dt, &newDT, sizeof(newDT));
	}
	float newK = singleTone->getFloat("K0", -1.0f);
	if (newK != -1.0f) {
		cudaMemcpyToSymbol(K0, &newK, sizeof(newK));
	}
}

void runFWMethod(simulationInput *simulation) {
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
		if(mkdirAndchdir("FW") != 1){
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
	curandInitialization << <simulation->blockSize, simulation->threadSize >> >(simulation->state);
	gpuErrchk(cudaDeviceSynchronize());
	int iterations = singleTone->getInt("bilions", 10);
	setConstants(singleTone);
	if (simulation->threadSize == 1024) {
		cudaFuncSetAttribute(trajectorySimulation, cudaFuncAttributeMaxDynamicSharedMemorySize, 65536);
	}
	for (int k = 0; k < iterations; ++k) {
		for (int i = 0; i < 31; ++i) {
			nullCount << <1, 1 >> >();
			gpuErrchk(cudaDeviceSynchronize());
			wCalc << <simulation->blockSize, simulation->threadSize >> >(simulation->w, simulation->pinj, i, simulation->history);
			gpuErrchk(cudaDeviceSynchronize());
			trajectorySimulation << <simulation->blockSize, simulation->threadSize, simulation->threadSize * sizeof(curandState_t) + 
				simulation->threadSize * sizeof(float2) >> >(simulation->pinj, simulation->history, i, simulation->state);
			gpuErrchk(cudaDeviceSynchronize());
			cudaMemcpyFromSymbol(&counter, countAtomic, sizeof(int), 0, cudaMemcpyDeviceToHost);
			if (counter != 0) {
				gpuErrchk(cudaMemcpy(simulation->local_history, simulation->history, counter*sizeof(trajectoryHistory), cudaMemcpyDeviceToHost));
				for (int j = 0; j < counter; ++j) {
					fmt::fprintf(bb, " %g  %g  %g  %g %g \n", simulation->pinj[simulation->local_history[j].id], simulation->local_history[j].p, simulation->local_history[j].r, simulation->w[simulation->local_history[j].id], simulation->local_history[j].sumac);
				}
			}
		}
	}
	fclose(bb);
}