#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <chrono>
#include <ctime>
#include <curand.h>
#include <curand_kernel.h>
#include <thrust/count.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/execution_policy.h>
#include <stack>

#include "../../Constants/CosmicConstants.cuh"
#include "../../Constants/ParamsCarrier.cuh"
#include "TwoDimensionFPDefinition.cuh"
#include "cudaErrorCheck.cuh"
#include "../../libs/FMT/format.h"
#include "../../libs/FMT/printf.h"

extern "C" void runTwoDimensionFPMethod(simulationInputTwoFP *simulation);
extern "C" int mkdirAndchdir(std::string directoryName);

__device__ int countAtomic;

__global__ void nullCountTwoFP() {
	countAtomic = 0;
}

__device__ float getTkininjTwoFP(unsigned long long state) {
	unsigned long long  ownState = state;
	int modulo;
	if (ownState > (151 * 10000)) {
		ownState -= (__double2ll_rd(ownState / (151 * 10000)) * (151 * 10000));
	}
	if (ownState >= 10000) {
		modulo = __double2int_rd(ownState / 10000);
	}
	else {
		modulo = 0;
	}
	return ((modulo)+((ownState - (modulo * 10000) + 1) / 10000.0f));
}

__global__ void curandInitializationTwoFP(curandState_t *state) {
	int execID = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(clock(), execID, 0, &state[execID]);
}

__global__ void wCalcTwoFP(float *p, float *mu, double *w, int padding, curandState* state) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int idx = threadIdx.x; 
	float Tkinw = getTkininjTwoFP(BLOCK_SIZE * THREAD_SIZE * padding + id)*1e9f*q;
	float Rig = sqrtf(Tkinw*(Tkinw + (2 * T0w))) / q;
    float pinj = Rig*q / c;
    p[id] = pinj; 
	double newW = (m0_double*m0_double*c_double*c_double*c_double*c_double) + (pinj*pinj*c_double*c_double);
	newW = (pow(newW, -1.85) / pinj) / 1e45;
	w[id] = newW;
	__shared__ curandState cuState[THREAD_SIZE];
	cuState[idx] = state[id];
	mu[id] = cosf(curand_uniform(&cuState[idx]) * Pi / 2); 
	state[id] = cuState[idx];
}

__global__ void trajectorySimulationTwoFP(float *pinj, float *muinj, trajectoryHistoryTwoFP* history, int padding, curandState* state) {
	extern __shared__ int sharedMemory[];
	int id = blockIdx.x * blockDim.x + threadIdx.x; 
	int idx = threadIdx.x;
	float r =  100.0001f;
    float mu = muinj[id]; 
    float p = pinj[id];
	float beta, Rig, dr, dp, tem1, Krr, dLBT, Kpar, dmu, dKrr, Ktt, Kper, oldMu, LBT = 0;
	float Tkin = getTkininjTwoFP(BLOCK_SIZE * THREAD_SIZE * padding + id);
	float2* generated = (float2*) sharedMemory;
	curandState* cuState = (curandState *)(&generated[THREAD_SIZE]);
	cuState[idx] = state[id];
	int count;
	for (; r < 100.0002f;) {
        beta = sqrtf(Tkin*(Tkin + T0 + T0)) / (Tkin + T0);
		Rig = sqrtf(Tkin*(Tkin + (2.0f * T0)));
        tem1 = (1.0f + (omega2*r*r*(1.0f-(mu * mu))/(V*V)));
        Kpar = K0*beta*Rig*(r*r / sqrtf(tem1)); 
        Kper = kparKper * Kpar;
        Krr = Kper + ((Kpar - (Kper))/tem1);
        Ktt = Kper;
		dKrr = -1.0f*((Kpar) - (Kper));
		dKrr *= (omega2*2.0f*r*(1.0f-(mu * mu))/(V*V))/(tem1 * tem1);
		generated[idx] = curand_box_muller(&cuState[idx]);
		dr = (V + (2*Krr/r) + dKrr)*dt + (generated[idx].x*sqrtf(2.0f*Krr*dt));
		dmu = mu*Ktt*dt/(r*r*sqrtf(1.0f-(mu*mu)));  
		dmu += generated[idx].y*sqrtf(2.0f*Ktt*dt/(r*r));
        
        dp = 2.0f*V*p*dt/(3.0f*r);
        dLBT = -4.0f*V*dt/(3.0f*r);

        r += dr;
        p -= dp;
		mu += dmu;
        oldMu = mu;
        LBT += dLBT;

        Rig = p*c / q;
		Tkin = (sqrtf((T0*T0*q*q*1e9f*1e9f) + (q*q*Rig *Rig)) - (T0*q*1e9f)) / (q*1e9f);
        beta = sqrtf(Tkin*(Tkin + T0 + T0)) / (Tkin + T0);
		
		if (mu>1.0f) {
            mu -= (2.0f*(fabs(mu) - 1.0f));
        } else if (mu<-1.0f) {
            mu += (2.0f*(fabs(mu) - 1.0f));
        }
            
		if (r<0.1f) {
            r -= dr;
            mu = oldMu; 
			p += dp;
		}
		if (beta>0.001f&&Tkin<200.0f) {
			if ((r - 1.0f) / ((r - dr) - 1.0f) < 0.0f) {
				count = atomicAdd(&countAtomic, 1);
				history[count].setValues(p, r, id, mu, LBT);
			}
		}
		else if (beta<0.001f) {
			break;
		}
	}
	state[id] = cuState[idx];
}


void setConstantsTwoFP(ParamsCarrier *singleTone) {
	float newDT = singleTone->getFloat("dt", -1.0f);
	if (newDT != -1.0f) {
		cudaMemcpyToSymbol(dt, &newDT, sizeof(newDT));
	}
	float newK = singleTone->getFloat("K0", -1.0f);
	if (newK != -1.0f) {
		cudaMemcpyToSymbol(K0, &newK, sizeof(newK));
	}
	float newKparKper = singleTone->getFloat("kparKper", -1.0f);
	if (newKparKper != -1.0f) {
		cudaMemcpyToSymbol(kparKper, &newKparKper, sizeof(newKparKper));
	}
}

void runTwoDimensionFPMethod(simulationInputTwoFP *simulation) {
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
		if(mkdirAndchdir("2DFP") != 1){
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
	curandInitializationTwoFP << <simulation->blockSize, simulation->threadSize >> >(simulation->state);
	gpuErrchk(cudaDeviceSynchronize());
	int iterations = singleTone->getInt("bilions", 10);
	setConstantsTwoFP(singleTone);
	if (simulation->threadSize == 1024) {
		cudaFuncSetAttribute(trajectorySimulationTwoFP, cudaFuncAttributeMaxDynamicSharedMemorySize, 65536);
	}
	#pragma unroll
	for (int k = 0; k < iterations; ++k) {
		#pragma unroll
		for (int i = 0; i < 31; ++i) {
			nullCountTwoFP << <1, 1 >> >();
			gpuErrchk(cudaDeviceSynchronize());
			wCalcTwoFP << <simulation->blockSize, simulation->threadSize >> >(simulation->pinj, simulation->mu, simulation->w, i, simulation->state);
			gpuErrchk(cudaDeviceSynchronize());
			trajectorySimulationTwoFP << <simulation->blockSize, simulation->threadSize >> >(simulation->pinj, simulation->mu, simulation->history, i, simulation->state);
			gpuErrchk(cudaDeviceSynchronize());
			cudaMemcpyFromSymbol(&counter, countAtomic, sizeof(int), 0, cudaMemcpyDeviceToHost);
			if (counter != 0) {
				gpuErrchk(cudaMemcpy(simulation->local_history, simulation->history, counter*sizeof(trajectoryHistoryTwoFP), cudaMemcpyDeviceToHost));
				for (int j = 0; j < counter ; ++j) {
					fmt::fprintf(bb, "%g %g %g %g %g %g\n", simulation->pinj[simulation->local_history[j].id] , simulation->local_history[j].p, simulation->local_history[j].r, simulation->local_history[j].mu, simulation->w[simulation->local_history[j].id], simulation->local_history[j].LBT);
				}				
			}
		}
	}
	fclose(bb);
}
