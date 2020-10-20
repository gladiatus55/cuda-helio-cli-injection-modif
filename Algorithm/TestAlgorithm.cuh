#ifndef TEST_ALGORITHM_H
#define TEST_ALGORITHM_H

#include "AbstractAlgorithm.cuh"

class TestAlgorithm : public AbstractAlgorithm {
public:
	void runAlgorithm(ParamsCarrier *singleTone);
};

#endif // !TEST_ALGORITHM_H
