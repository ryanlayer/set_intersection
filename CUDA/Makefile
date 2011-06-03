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
LFLAGS_A=-fPIC
LFLAGS_B=-L/usr/local/cuda/lib64 \
		-L$(SDK_PATH)/C/lib \
		-L$(SDK_PATH)/C/common/lib/linux \
		-L$(SDK_PATH)/shared/lib \
		-lcudart \
		-lcudpp_x86_64 

NVCC_OPTS=-DUNIX -O2

SET_FILE_OBJS=../lib/set_intersect.o
TIMER_FILE_OBJS=../lib/timer.o
BED_FILE_OBJS=../lib/bed.o

OBJS= $(SET_FILE_OBJS) $(TIMER_FILE_OBJS) $(BED_FILE_OBJS)

SET_FILE_SRC=../lib/set_intersect.c
TIMER_FILE_SRC=../lib/timer.c
BED_FILE_SRC=../lib/bed.c
MT_SRC=../lib/mt.c

SRC= $(SET_FILE_SRC) $(TIMER_FILE_SRC) $(BED_FILE_SRC) $(MT_SRC)

EXECS=cuda_bsearch \
	  test_cuda \
	  num_sim_bsearch_gm_cuda \
	  num_sim_bsearch_sm_cuda \
	  count_bsearch_gm_cuda \
	  big_count_bsearch_gm_cuda \
	  enumerate_bsearch_gm_cuda \
	  count_bsearch_sm_cuda \
	  num_sim_brute_force_gm_cuda\
	  num_sim_test \
	  count_brute_force_gm_cuda \
	  bsearch_gm_cuda \
	  pre_gm_index_wbsearch_gm_cuda

.SUFFIXES: .c .cpp .cu .o

cuda_bsearch: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  cuda_bsearch.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o\
		cuda_bsearch.cu.o \
		$(LFLAGS_B)

bsearch_gm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  bsearch_gm_cuda.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o\
		bsearch_gm_cuda.cu.o \
		$(LFLAGS_B)

pre_gm_index_wbsearch_gm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  pre_gm_index_wbsearch_gm_cuda.cu.o
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o\
		pre_gm_index_wbsearch_gm_cuda.cu.o \
		$(LFLAGS_B)



num_sim_bsearch_gm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  num_sim_bsearch_gm_cuda.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o \
		num_sim_bsearch_gm_cuda.cu.o \
		$(LFLAGS_B)


num_sim_test: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  num_sim_test.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o \
	    num_sim_test.cu.o \
		$(LFLAGS_B)

num_sim_brute_force_gm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  num_sim_brute_force_gm_cuda.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o \
	  	num_sim_brute_force_gm_cuda.cu.o \
		$(LFLAGS_B)

count_brute_force_gm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  count_brute_force_gm_cuda.cu.o
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o \
		count_brute_force_gm_cuda.cu.o \
		$(LFLAGS_B)


count_bsearch_gm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  count_bsearch_gm_cuda.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o \
		count_bsearch_gm_cuda.cu.o \
		$(LFLAGS_B)

big_count_bsearch_gm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  big_count_bsearch_gm_cuda.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o \
		big_count_bsearch_gm_cuda.cu.o \
		$(LFLAGS_B)


count_bsearch_sm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  count_bsearch_sm_cuda.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o \
		count_bsearch_sm_cuda.cu.o \
		$(LFLAGS_B)

enumerate_bsearch_gm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  enumerate_bsearch_gm_cuda.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o \
		enumerate_bsearch_gm_cuda.cu.o \
		$(LFLAGS_B)

num_sim_bsearch_sm_cuda: radixsort.cu.o \
			  radixsort.o \
			  random.cu.o \
			  set_intersect_cuda.cu.o \
			  num_sim_bsearch_sm_cuda.cu.o 
	$(CPP) $(LFLAGS_A) -o $@ \
		$(SRC) \
		radixsort.cu.o \
		random.cu.o \
		radixsort.o \
		set_intersect_cuda.cu.o \
		num_sim_bsearch_sm_cuda.cu.o \
		$(LFLAGS_B)

test_cuda: set_intersect_cuda.cu.o \
		   test_cuda.cu.o
	$(CPP) $(LFLAGS_A) -o $@ \
		set_intersect_cuda.cu.o \
		test_cuda.cu.o \
		$(LFLAGS_B)

%.cu.o: %.cu
	$(NVCC)  $(INCLUDE) -c $< -o $@

%.o: %.cpp
	$(CPP)  $(INCLUDE) -c $< -o $@

clean:
	rm -f $(EXECS) *.a *.o $(LIBS) $(OBJS)
