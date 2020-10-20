#include "CosmicFactory.cuh"

#include "../Algorithm/TestAlgorithm.cuh"
#include "../Algorithm/FloatFWMethod.cuh"
#include "../Algorithm/FloatBPMethod.cuh"
#include "../Algorithm/TwoDimensionFPMethod.cuh"

AbstractAlgorithm* CosmicFactory::getAlgorithm(std::string name) {
	if (name.compare("TestRun") == 0) {
		return new TestAlgorithm();
	}else if (name.compare("FWMethod") == 0) {
		return new FloatFWMethod(); 
	}else if (name.compare("BPMethod") == 0) {
		return new FloatBPMethod(); 
	}else if (name.compare("TwoDimensionFPMethod") == 0) {
		return new TwoDimensionFPMethod();
	}
	return NULL; 
}