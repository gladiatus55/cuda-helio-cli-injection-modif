#ifndef ABSTRACT_ALGORITHM_H
#define ABSTRACT_ALGORITHM_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>
#include "../Constants/ParamsCarrier.cuh"

class AbstractAlgorithm {
public:
	virtual void runAlgorithm(ParamsCarrier*singleTone) = 0;
};

#endif // !ABSTRACT_ALGORITHM_H
