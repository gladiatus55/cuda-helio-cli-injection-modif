#ifndef FLOAT_FW_RESULT_H
#define FLOAT_FW_RESULT_H

#include "AbstractAlgorithm.cuh"

class FWMethodResult : public AbstractAlgorithm {
public:
	void runAlgorithm(ParamsCarrier *singleTone);
private: 
	int charcount( FILE *const fin ); 
};

#endif
