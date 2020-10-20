#ifndef COSMIC_FACTORY_H
#define COSMIC_FACTORY_H

#include "AbstractAlgorithmFactory.cuh"

class CosmicFactory : public AbstractAlgorithmFactory {
public: 
	AbstractAlgorithm* getAlgorithm(std::string name); 
};

#endif // !COSMIC_FACTORY_H
