#ifndef PARSE_PARAMS_H
#define PARSE_PARAMS_H

#include <string>
#include "ParamsCarrier.cuh"

class ParseParams{
public:
	int parseParams(int argc, char **argv);
	ParamsCarrier *getParams(); 
private: 
	ParamsCarrier *singleTone;
};

#endif