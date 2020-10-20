#ifndef FLOAT_TWO_FP_H
#define FLOAT_TWO_FP_H

#include <stdio.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include "AbstractAlgorithm.cuh"
#include "TwoDimensionFPResult.cuh"
#include "cuda_implementation/TwoDimensionFPDefinition.cuh"

extern "C" void runTwoDimensionFPMethod(simulationInputTwoFP *simulation);

class TwoDimensionFPMethod : public AbstractAlgorithm {
public:
	void runAlgorithm(ParamsCarrier *singleTone);
private: 
	int threadSize, blockSize; 
	void setThreadBlockSize(); 
};

#endif