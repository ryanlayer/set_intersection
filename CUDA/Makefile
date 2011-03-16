CPP=g++
SDK_PATH=$(HOME)/NVIDIA_GPU_Computing_SDK
SDK_COMMON_INC=$(SDK_PATH)/C/common/inc
SDK_SHARED_INC=$(SDK_PATH)/C/shared/inc
CUDA_INCLUDE=/usr/local/cuda/include
NVCC=/usr/local/cuda/bin/nvcc
INCLUDE=-I. \
		-I$(CUDA_INCLUDE) \
		-I$(SDK_COMMON_INC) \
	   	-I$(SDK_SHARED_INC)
CFLAGS=
LFLAGS=-fPIC -L/usr/local/cuda/lib64 \
		-L$(SDK_PATH)/C/lib \
		-L$(SDK_PATH)/C/common/lib/linux \
		-L$(SDK_PATH)/shared/lib \
		-lcudart \
		-lcuda \
		-lcudpp_x86_64

NVCC_OPTS=-DUNIX -O2

SET_FILE_OBJS=../lib/set_intersect.o
TIMER_FILE_OBJS=../lib/timer.o
BED_FILE_OBJS=../lib/bed.o

SET_FILE_SRC=../lib/set_intersect.c
TIMER_FILE_SRC=../lib/timer.c
BED_FILE_SRC=../lib/bed.c

EXECS=cuda_bsearch

OBJS= $(SET_FILE_OBJS) $(TIMER_FILE_OBJS) $(BED_FILE_OBJS)
SRC= $(SET_FILE_SRC) $(TIMER_FILE_SRC) $(BED_FILE_SRC)

.SUFFIXES: .c .cpp .cu .o

#cuda_bsearch: radixsort.cu.o random.cu.o cuda_bsearch.cu.o
cuda_bsearch: radixsort.cu.o random.cu.o cuda_bsearch.cu.o radixsort.o
	$(CPP) -fPIC -o $@ \
		$(SRC) \
		radixsort.cu.o random.cu.o cuda_bsearch.cu.o radixsort.o\
		-L/usr/local/cuda/lib64 \
		-L$(SDK_PATH)/C/lib \
		-L$(SDK_PATH)/C/common/lib/linux \
		-L$(SDK_PATH)//shared/lib \
		-lcudart \
		-lcudpp_x86_64 \
		-L/usr/local/cuda/lib64 




%.cu.o: %.cu
	$(NVCC)  $(INCLUDE) -c $< -o $@

%.o: %.cpp
	$(CPP)  $(INCLUDE) -c $< -o $@

clean:
	rm -f $(EXECS) *.a *.o $(LIBS) $(OBJS)

#g++ -fPIC   -o ../../bin/linux/release/cuda_bsearch obj/x86_64/release/radixsort.cpp.o bed.c set_intersect.c timer.c  obj/x86_64/release/cuda_bsearch.cu.o obj/x86_64/release/radixsort.cu.o obj/x86_64/release/random.cu.o      -L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart    -lcudpp_x86_64 -L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart -lcutil_x86_64 -lshrutil_x86_64


