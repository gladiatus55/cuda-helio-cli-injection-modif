#include "TestAlgorithm.cuh"
#include "../Constants/CosmicConstants.cuh"
#include "../Constants/ParamsCarrier.cuh"
#include "cuda_implementation/cudaErrorCheck.cuh"

void TestAlgorithm::runAlgorithm(ParamsCarrier *singleTone) {
	printf("I am alive!\n");
	float newDT = singleTone->getFloat("dt", -1.0f);
	cudaMemcpyToSymbol(dt, &newDT, sizeof(newDT));
	gpuErrchk(cudaPeekAtLastError());
	printf("New value of dt is %f!\n", newDT);
	printf("Everything is alright!\n");
}