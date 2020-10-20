#ifndef CUDA_ERROR_CHECK_H
#define CUDA_ERROR_CHECK_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPU fail: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) {
			exit(code);
		}
	}
}

#endif // !CUDA_ERROR_CHECK_H
