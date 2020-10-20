#ifndef BP_METHOD_H
#define BP_METHOD_H

#include "AbstractAlgorithm.cuh"
#include "cuda_implementation/BPDefines.cuh"
#include <stdio.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>

extern "C" void runBPMethod(simulationInputBP *simulation);

class FloatBPMethod : public AbstractAlgorithm {
public:
	void runAlgorithm(ParamsCarrier *singleTone);
private: 
	int threadSize, blockSize; 
	void setThreadBlockSize(); 
};

#endif




