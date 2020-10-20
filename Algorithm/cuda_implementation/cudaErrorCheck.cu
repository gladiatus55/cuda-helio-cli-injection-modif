#include "cudaErrorCheck.cuh"
#include <string>

extern "C" int mkdirAndchdir(std::string directoryName);

int mkdirAndchdir(std::string directoryName){
	char *name = new char[directoryName.length() + 1]; 
	strcpy(name , directoryName.c_str());
	DIR* dir = opendir(name);
	if (dir){
		int resultMkdir = chdir(name); 
		if (resultMkdir  != 0) {
    			printf("ERROR %d: unable to change directory; %s\n", resultMkdir , strerror(resultMkdir ));
			delete[] name;
			return -1; 
  		}	
	}else if (ENOENT == errno){		
		int resultMkdir = mkdir(name,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (resultMkdir  != 0) {
    			printf("ERROR %d: unable to mkdir; %s\n", resultMkdir , strerror(resultMkdir ));
			delete[] name;
			return -1; 
  		}
		resultMkdir = chdir(name); 
		if (resultMkdir  != 0) {
    			printf("ERROR %d: unable to change directory; %s\n", resultMkdir , strerror(resultMkdir ));
			delete[] name;
			return -1; 
  		}		
	}else{
		printf("ERROR: unable to open directory");
		delete[] name;
		return -1; 
	}	
	delete[] name;
	return 1; 
}
