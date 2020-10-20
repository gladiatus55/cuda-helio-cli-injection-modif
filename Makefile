NVCC = nvcc
ifeq ($(turing),1)
    GENCODE = -gencode=arch=compute_75,code=\"sm_75,compute_75\"
else
    GENCODE = -gencode=arch=compute_61,code=\"sm_61,compute_61\"
endif
CFLAGS =--compiler-options -w -w -O3 -Xptxas -O3 -Xcompiler -O3 --use_fast_math -std=c++11

csrc = $(wildcard Constants/*.cu) \
       $(wildcard Algorithm/*.cu) \
       $(wildcard Algorithm/cuda_implementation/*.cu) \
       $(wildcard Factories/*.cu) \
       $(wildcard *.cu) 
ccsrc = $(wildcard libs/FMT/*.cc) 
obj = $(csrc:.cu=.o) $(ccsrc:.cc=.o)

cosmicCUDA: $(obj)
	$(NVCC) $(GENCODE) -o $@ $^ $(CFLAGS)

clean:
	rm -f $(obj) all

%.o: %.cu
	$(NVCC) $(GENCODE) $(CFLAGS) -o $@ -c $<