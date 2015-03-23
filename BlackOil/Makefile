INC_ADETL= -I../../../Library/ADETL
#INC_SUPERLU= -I/opt/SuperLU_SERIAL/include/
#INC= -I./src $(INC_ADETL) $(INC_SUPERLU)
#MKLRTDIR = /opt/intel/mkl/lib/intel64
#LIBS= -lsuperlu -L$(MKLRTDIR) -Wl,--start-group $(MKLRTDIR)/libmkl_intel_lp64.a $(MKLRTDIR)/libmkl_core.a $(MKLRTDIR)/libmkl_sequential.a -Wl,--end-group -lpthread -lm
#LIBDIR= -L/opt/SuperLU_SERIAL/lib
CFLAGS= -Wall -g -O3 -std=c++11
SRC= *.cpp 

all:
	g++ $(INC_ADETL) $(CFLAGS) $(SRC) -o ./run

clean:
	rm -rf ./run ./*~
