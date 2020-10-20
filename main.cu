#include <stdio.h>
#include <iostream>

#include "Constants/CosmicConstants.cuh"
#include "Constants/ParamsCarrier.cuh"
#include "Constants/ParseParams.cuh"
#include "Algorithm/AbstractAlgorithm.cuh"
#include "Factories/AbstractAlgorithmFactory.cuh"

int main(int argc, char **argv){
	AbstractAlgorithmFactory *factory = AbstractAlgorithmFactory::CreateFactory
		(AbstractAlgorithmFactory::TYPE_ALGORITHM::COSMIC);
	ParseParams *parse = new ParseParams(); 
	if(parse->parseParams(argc, argv) != 1){
		return -1; 
	} 		
	ParamsCarrier*singleTone;
	singleTone = parse->getParams();
	AbstractAlgorithm *actualAlgorithm;
	actualAlgorithm = factory->getAlgorithm(singleTone->getString("algorithm","FWMethod"));
	actualAlgorithm->runAlgorithm(singleTone);
	return 0;
}