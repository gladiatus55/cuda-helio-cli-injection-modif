#ifndef BP_RES_H
#define BP_RES_H

#include "AbstractAlgorithm.cuh"

class BPResult : public AbstractAlgorithm {
public:
	void runAlgorithm(ParamsCarrier *singleTone);
private: 
	int charcount( FILE *const fin ); 
};

#endif
