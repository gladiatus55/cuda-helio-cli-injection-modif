#ifndef FLOAT_FW_H
#define FLOAT_FW_H

#include <stdio.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include "AbstractAlgorithm.cuh"
#include "FWMethodResult.cuh"
#include "cuda_implementation/FloatFWDefines.cuh"

extern "C" void runFWMethod(simulationInput *simulation);

class FloatFWMethod : public AbstractAlgorithm {
public:
	void runAlgorithm(ParamsCarrier *singleTone);
private: 
	int threadSize, blockSize; 
	void setThreadBlockSize(); 
};

#endif


