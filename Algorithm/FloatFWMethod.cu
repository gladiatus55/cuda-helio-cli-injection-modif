#include "FloatFWMethod.cuh"
#include "FWMethodResult.cuh"
#include "cuda_implementation/cudaErrorCheck.cuh"

void FloatFWMethod::runAlgorithm(ParamsCarrier *singleTone){
	simulationInput simulation; 
	setThreadBlockSize(); 
	curandState_t *state;
	double *w;
	float *pinj;
	trajectoryHistory *history, *local_history;
	gpuErrchk( cudaMallocManaged(&w, ((blockSize * threadSize) *sizeof(double))) );
	gpuErrchk( cudaMallocManaged(&pinj, ((blockSize * threadSize) *sizeof(float))) );
	gpuErrchk( cudaMallocManaged(&state, ((blockSize * threadSize) * sizeof(curandState_t))) );
	gpuErrchk(cudaMallocHost(&local_history, ((BLOCK_SIZE * THREAD_SIZE * 10) *sizeof(trajectoryHistory))));
	gpuErrchk(cudaMalloc(&history, ((BLOCK_SIZE * THREAD_SIZE * 10) *sizeof(trajectoryHistory))));
	simulation.singleTone = singleTone;
	simulation.history = history; 
	simulation.local_history = local_history; 
	simulation.pinj = pinj; 
	simulation.state = state; 
	simulation.w = w; 
	simulation.threadSize = threadSize; 
	simulation.blockSize = blockSize; 
	runFWMethod(&simulation);
	gpuErrchk( cudaFree(w) );
	gpuErrchk( cudaFree(pinj) );
	gpuErrchk( cudaFree(state) );
	gpuErrchk( cudaFree(history) );
	gpuErrchk( cudaFreeHost(local_history) );
	AbstractAlgorithm *result;
	result = new FWMethodResult(); 
	result->runAlgorithm(singleTone);
}

// Compute capability actual device 
void FloatFWMethod::setThreadBlockSize() {
	cudaDeviceProp gpuProperties; 
	gpuErrchk( cudaGetDeviceProperties(&gpuProperties, 0) );
	int computeCapability = gpuProperties.major * 100 + gpuProperties.minor * 10; 
	switch (computeCapability) {
		case 610:	blockSize = 65536;
				threadSize = 512;
				break;
		case 750:	blockSize = 32768;
				threadSize = 1024;
				break;
		default:	blockSize = 64;
				threadSize = 64;
				break;
	}
}