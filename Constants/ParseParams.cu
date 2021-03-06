#include "ParseParams.cuh"
#include <string>
#include "../libs/CLI/CLI.hpp"
#include "../libs/CLI/App.hpp"
#include "../libs/CLI/Option.hpp"

#include "CosmicConstants.cuh"
#include "ParamsCarrier.cuh"

int ParseParams::parseParams(int argc, char **argv){
	float newDT, newK, newV, newKparKper;
	std::string newDestination;
	CLI::Option* backwardMethod,*runTest, *forwardMethod, *csv, *twoDimensionFP;
	int bilions; 
	singleTone = singleTone->instance();
	CLI::App app{ "App description" };
	runTest = app.add_flag("-t,--test", "Run a short test");
	forwardMethod = app.add_flag("-F,--Forward", "Run a forward method(base - 1 bilion)")->group("Methods");
	backwardMethod = app.add_flag("-B,--Backward", "Run a backward method(base - 16 milions)")->group("Methods");	
	twoDimensionFP = app.add_flag("-P,--TwoFP", "Run a 2D Fp method(base - 1 bilion)")->group("Methods");	
	csv = app.add_flag("-c,--csv", "Output will be in .csv");
	CLI::Option* dtset = app.add_option("-d,--dt", newDT, "Set dt to new value");
	CLI::Option* kset = app.add_option("-K,--K0", newK ,"Set K to new value(cm^2/s)");
	CLI::Option* vset = app.add_option("-V,--V", newV ,"Set V to new value(km/s)");
	CLI::Option* kparKper = app.add_option("-E,--KparKper", newKparKper ,"Set KparKper");
	CLI::Option* destination = app.add_option("-p,--path", newDestination ,"Set destination folder name");
	CLI::Option* setBilions = app.add_option("-N,--Bilions", bilions, "Set number of iteration in base units for each method");
	backwardMethod->excludes(forwardMethod); 
	backwardMethod->excludes(runTest); 
	backwardMethod->excludes(twoDimensionFP); 
	forwardMethod->excludes(backwardMethod);
	forwardMethod->excludes(twoDimensionFP); 
	forwardMethod->excludes(runTest);
	runTest->excludes(backwardMethod); 
	runTest->excludes(twoDimensionFP); 
	runTest->excludes(forwardMethod); 
	runTest->excludes(csv); 
	twoDimensionFP->excludes(backwardMethod);
	twoDimensionFP->excludes(forwardMethod);
	twoDimensionFP->excludes(runTest);
	CLI11_PARSE(app, argc, argv);
	if(!*forwardMethod && !*backwardMethod && !*runTest && !*twoDimensionFP){			
		printf("At least one method must be selected!"); 
		return -1; 
	}
	if(*csv){
		singleTone->putInt("csv", 1);
	}
	if (*dtset) {
		if(newDT< 0.1f || newDT> 1000){
			printf("dt is out of range!(0.1-1000)"); 
			return -1; 
		}
		singleTone->putFloat("dt",newDT);
	}
	if (*kset) {
		if(newK < 1e22 || newK > 1e23){
			printf("K0 is out of range!(1e22-1e23 cm^2/s)"); 
			return -1; 
		}
		char buffer[80];
   		sprintf(buffer, "%g", newK);
		singleTone->putString("K0input", buffer); 
		newK = newK * 4.4683705e-27;
		singleTone->putFloat("K0", newK); 
	}
	if (*setBilions) {
		if(bilions <= 0){
			printf("Number of iteration must be greater than 0!"); 
			return -1; 
		}
		singleTone->putInt("bilions", bilions);
	}
	if (*kparKper) {
		if(newKparKper > 0.1 && newKparKper < 0.01){
			printf("KparKper is out of range(0.01-0.1)!"); 
			return -1; 
		}
		singleTone->putInt("kparKper", newKparKper);
	}
	if(*vset){
		if(newV < 100 || newV > 1500){
			printf("V is out of range!(100-1500 km/s)"); 
			return -1; 
		}
		singleTone->putString("Vinput", std::to_string(newV)); 
		newV = newV * 6.68458712e-9; 
		singleTone->putFloat("V", newV); 
	}
	if(*destination) {
		singleTone->putString("destination", newDestination);
	}
	if (*forwardMethod) {
		singleTone->putString("algorithm","FWMethod");
	}else if(*backwardMethod){
		singleTone->putString("algorithm","BPMethod");	
	}else if (*runTest) {
		singleTone->putString("algorithm","TestRun");
	}else if (*twoDimensionFP) {
		singleTone->putString("algorithm","TwoDimensionFPMethod");
	}
	return 1; 
}

ParamsCarrier* ParseParams::getParams(){
	return singleTone; 
}
