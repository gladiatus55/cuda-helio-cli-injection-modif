#ifndef ABSTRACT_ALGORITHM_FACTORY_H
#define ABSTRACT_ALGORITHM_FACTORY_H

#include <string>

#include "../Algorithm/AbstractAlgorithm.cuh"

class AbstractAlgorithmFactory {
public: 
	enum TYPE_ALGORITHM{
		COSMIC
	};

	virtual AbstractAlgorithm* getAlgorithm(std::string name) = 0; 

	static AbstractAlgorithmFactory* CreateFactory(TYPE_ALGORITHM factory);
};

#endif // !ABSTRACT_ALGORITHM_FACTORY
