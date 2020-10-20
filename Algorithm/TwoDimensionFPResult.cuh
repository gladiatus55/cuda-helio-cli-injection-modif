#ifndef FLOAT_TWO_FP_RESULT_H
#define FLOAT_TWO_FP_RESULT_H

#include "AbstractAlgorithm.cuh"

class TwoDimensionFPResult : public AbstractAlgorithm {
public:
	void runAlgorithm(ParamsCarrier *singleTone);
private: 
	int charcount( FILE *const fin ); 
};

#endif
