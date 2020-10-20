#include "TwoDimensionFPMethod.cuh"
#include "cuda_implementation/cudaErrorCheck.cuh"

void TwoDimensionFPMethod::runAlgorithm(ParamsCarrier *singleTone){
	simulationInputTwoFP simulation; 
	setThreadBlockSize(); 
	curandState_t *state;
	double *w;
    float *pinj, *mu;
	trajectoryHistoryTwoFP *history, *local_history;
	gpuErrchk( cudaMallocManaged(&w, ((blockSize * threadSize) *sizeof(double))) );
	gpuErrchk( cudaMallocManaged(&mu, ((blockSize * threadSize) *sizeof(float))) );
	gpuErrchk( cudaMallocManaged(&pinj, ((blockSize * threadSize) *sizeof(float))) );
	gpuErrchk( cudaMallocManaged(&state, ((blockSize * threadSize) * sizeof(curandState_t))) );
	gpuErrchk(cudaMallocHost(&local_history, ((BLOCK_SIZE * THREAD_SIZE * 10) *sizeof(trajectoryHistoryTwoFP))));
	gpuErrchk(cudaMalloc(&history, ((BLOCK_SIZE * THREAD_SIZE * 10) *sizeof(trajectoryHistoryTwoFP))));
	simulation.singleTone = singleTone;
	simulation.history = history; 
	simulation.local_history = local_history; 
	simulation.pinj = pinj; 
	simulation.state = state; 
    simulation.w = w; 
    simulation.mu = mu;
	simulation.threadSize = threadSize; 
	simulation.blockSize = blockSize; 
	runTwoDimensionFPMethod(&simulation);
	gpuErrchk(cudaFree(w));
	gpuErrchk(cudaFree(pinj));
    gpuErrchk(cudaFree(state));
    gpuErrchk(cudaFree(mu))
	gpuErrchk(cudaFree(history));
	gpuErrchk(cudaFreeHost(local_history));
	AbstractAlgorithm *result;
	result = new TwoDimensionFPResult(); 
	result->runAlgorithm(singleTone);
}

// Compute capability actual device 
void TwoDimensionFPMethod::setThreadBlockSize() {
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